!> \file
!> $Id: opencmiss.f90 542 2009-06-03 17:16:22Z chrispbradley $
!> \author Chris Bradley
!> \brief The top level OpenCMISS module.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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
!>
!> The top level OpenCMISS module. This module is the buffer module between the OpenCMISS library and user code.
MODULE OPENCMISS

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CMISS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES
   
  IMPLICIT NONE

  PRIVATE
  
  !Module parameters
  
  !Module types

  !>Contains information about a basis function.
  TYPE CMISSBasisType
    PRIVATE
    TYPE(BASIS_TYPE), POINTER :: BASIS
  END TYPE CMISSBasisType

  !>Contains information on the boundary conditions for the equations set.
  TYPE CMISSBoundaryConditionsType
    PRIVATE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
  END TYPE CMISSBoundaryConditionsType

  !>Contains information on a control loop.
  TYPE CMISSControlLoopType
    PRIVATE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
  END TYPE CMISSControlLoopType

  !>Contains information on a coordinate system.
  TYPE CMISSCoordinateSystemType
    PRIVATE
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  END TYPE CMISSCoordinateSystemType

  !>Contains information on the mesh decomposition.
  TYPE CMISSDecompositionType
    PRIVATE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  END TYPE CMISSDecompositionType
  
  !>Contains information about the equations in an equations set.
  TYPE CMISSEquationsType
    PRIVATE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  END TYPE CMISSEquationsType
  
  !>Contains information on an equations set defined on a region. 
  TYPE CMISSEquationsSetType
    PRIVATE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  END TYPE CMISSEquationsSetType

  !>Contains information for a field defined on a region.
  TYPE CMISSFieldType
    PRIVATE
    TYPE(FIELD_TYPE), POINTER :: FIELD
  END TYPE CMISSFieldType
  
  !>Contains information for a fields defined on a region.
  TYPE CMISSFieldsType
    PRIVATE
    TYPE(FIELDS_TYPE), POINTER :: FIELDS
  END TYPE CMISSFieldsType
  
  !>Contains information on a generated mesh.
  TYPE CMISSGeneratedMeshType
    PRIVATE
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
  END TYPE CMISSGeneratedMeshType
  
  !>Contains information about a history file for a control loop.
  TYPE CMISSHistoryType
    PRIVATE
    TYPE(HISTORY_TYPE), POINTER :: HISTORY
  END TYPE CMISSHistoryType
  
  !>Contains information on a mesh defined on a region.
  TYPE CMISSMeshType
    PRIVATE
    TYPE(MESH_TYPE), POINTER :: MESH
  END TYPE CMISSMeshType
  
  !>Contains information on the nodes defined on a region.
  TYPE CMISSNodesType
    PRIVATE
    TYPE(NODES_TYPE), POINTER :: NODES
  END TYPE CMISSNodesType
  
  !>Contains information for a problem.
  TYPE CMISSProblemType
    PRIVATE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  END TYPE CMISSProblemType

  !>Contains information for a particular quadrature scheme for a basis. 
  TYPE CMISSQuadratureType
    PRIVATE
    TYPE(QUADRATURE_TYPE), POINTER :: QUADRATURE
  END TYPE CMISSQuadratureType

 !>Contains information for a region.
  TYPE CMISSRegionType
    PRIVATE
    TYPE(REGION_TYPE), POINTER :: REGION
  END TYPE CMISSRegionType

  !>Contains information about a solver.
  TYPE CMISSSolverType
    PRIVATE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
  END TYPE CMISSSolverType
  
  !>Contains information about the solver equations for a solver.
  TYPE CMISSSolverEquationsType
    PRIVATE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  END TYPE CMISSSolverEquationsType
  
  !Module variables

  TYPE(VARYING_STRING) :: ERROR

  !Interfaces

  INTERFACE CMISSInitialise
    MODULE PROCEDURE CMISSInitialiseNumber
    MODULE PROCEDURE CMISSInitialiseObj
  END INTERFACE !CMISSInitialise

  PUBLIC CMISSFinalise,CMISSInitialise

  PUBLIC CMISSBasisType,CMISSBasisTypeFinalise,CMISSBasisTypeInitialise

  PUBLIC CMISSBoundaryConditionsType,CMISSBoundaryConditionsTypeFinalise,CMISSBoundaryConditionsTypeInitialise

  PUBLIC CMISSControlLoopType,CMISSControlLoopTypeFinalise,CMISSControlLoopTypeInitialise

  PUBLIC CMISSCoordinateSystemType,CMISSCoordinateSystemTypeFinalise,CMISSCoordinateSystemTypeInitialise

  PUBLIC CMISSDecompositionType,CMISSDecompositionTypeFinalise,CMISSDecompositionTypeInitialise

  PUBLIC CMISSEquationsType,CMISSEquationsTypeFinalise,CMISSEquationsTypeInitialise

  PUBLIC CMISSEquationsSetType,CMISSEquationsSetTypeFinalise,CMISSEquationsSetTypeInitialise

  PUBLIC CMISSFieldType,CMISSFieldTypeFinalise,CMISSFieldTypeInitialise

  PUBLIC CMISSFieldsType,CMISSFieldsTypeFinalise,CMISSFieldsTypeInitialise

  PUBLIC CMISSGeneratedMeshType,CMISSGeneratedMeshTypeFinalise,CMISSGeneratedMeshTypeInitialise

  PUBLIC CMISSHistoryType,CMISSHistoryTypeFinalise,CMISSHistoryTypeInitialise

  PUBLIC CMISSMeshType,CMISSMeshTypeFinalise,CMISSMeshTypeInitialise

  PUBLIC CMISSNodesType,CMISSNodesTypeFinalise,CMISSNodesTypeInitialise

  PUBLIC CMISSProblemType,CMISSProblemTypeFinalise,CMISSProblemTypeInitialise

  PUBLIC CMISSQuadratureType,CMISSQuadratureTypeFinalise,CMISSQuadratureTypeInitialise

  PUBLIC CMISSRegionType,CMISSRegionTypeFinalise,CMISSRegionTypeInitialise

  PUBLIC CMISSSolverType,CMISSSolverTypeFinalise,CMISSSolverTypeInitialise

  PUBLIC CMISSSolverEquationsType,CMISSSolverEquationsTypeFinalise,CMISSSolverEquationsTypeInitialise
  

!!==================================================================================================================================
!!
!! ANALYTIC_ANALYSIS_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  !>Output the analytic error analysis for a field compared to the analytic values parameter set.
  INTERFACE CMISSAnalyticAnalysisOutput
    MODULE PROCEDURE CMISSAnalyticAnalysisOutputNumber
    MODULE PROCEDURE CMISSAnalyticAnalysisOutputObj
  END INTERFACE !CMISSAnalyticAnalysisOutput
  
  PUBLIC CMISSAnalyticAnalysisOutput
 
!!==================================================================================================================================
!!
!! BASE_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !> \addtogroup OPENCMISS_DiagnosticAndTimingConstants OPENCMISS::DiagnosticAndTiming::Constants
  !> \brief Diagnostic and Timing constants.
  !>@{  
  !> \addtogroup OPENCMISS_DiagnosticTypes OPENCMISS::DiagnosticTypes
  !> \brief Diganostic constants.
  !> \see OPENCMISS::DiagnosticAndTiming,OPENCMISS
  !>@{  
  INTEGER(INTG), PARAMETER :: CMISSAllDiagType = ALL_DIAG_TYPE !<Type for setting diagnostic output in all routines \see OPENCMISS_DiagnosticTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSInDiagType = IN_DIAG_TYPE !<Type for setting diagnostic output in one routine \see OPENCMISS_DiagnosticTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFromDiagType = FROM_DIAG_TYPE !<Type for setting diagnostic output in one routine downwards \see OPENCMISS_DiagnosticTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_TimingTypes OPENCMISS::TimingTypes
  !> \brief Timing constants.
  !> \see OPENCMISS::DiagnosticAndTiming,OPENCMISS
  !>@{  
  INTEGER(INTG), PARAMETER :: CMISSAllTimingType = ALL_TIMING_TYPE !<Type for setting timing output in all routines \see OPENCMISS_TimingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSInTimingType = IN_TIMING_TYPE !<Type for setting timing output in one routine \see OPENCMISS_TimingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFromTimingType = FROM_TIMING_TYPE !<Type for setting timing output from one routine downwards \see OPENCMISS_TimingTypes,OPENCMISS
  !>@}
  !>@}
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSAllDiagType,CMISSInDiagType,CMISSFromDiagType

  PUBLIC CMISSAllTimingType,CMISSInTimingType,CMISSFromTimingType

  PUBLIC CMISSDiagnosticsSetOff,CMISSDiagnosticsSetOn

  PUBLIC CMISSOutputSetOff,CMISSOutputSetOn

  PUBLIC CMISSTimingSetOff,CMISSTimingSetOn,CMISSTimingSummaryOutput

!!==================================================================================================================================
!!
!! BASIS_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !> \addtogroup OPENCMISS_BasisConstants OPENCMISS::BasisConstants
  !> \brief Basis function constants.
  !>@{  
  !> \addtogroup OPENCMISS_BasisTypes OPENCMISS::Basis::BasisTypes
  !> \brief Basis definition type parameters.
  !> \see OPENCMISS::BasisConstants,OPENCMISS
  !>@{ 
  INTEGER(INTG), PARAMETER :: CMISSBasisLagrangeHermiteTPType = BASIS_LAGRANGE_HERMITE_TP_TYPE !<Lagrange-Hermite tensor product basis type \see OPENCMISS_BasisTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisSimplexType = BASIS_SIMPLEX_TYPE !<Simplex basis type \see OPENCMISS_BasisTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisSerendipityType = BASIS_SERENDIPITY_TYPE !<Serendipity basis type \see OPENCMISS_BasisTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisAuxilliaryType = BASIS_AUXILLIARY_TYPE !<Auxillary basis type \see OPENCMISS_BasisTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisBSplineTPType = BASIS_B_SPLINE_TP_TYPE !<B-spline basis type \see OPENCMISS_BasisTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisFourierLagrangeHermiteTPType = BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE !<Fourier-Lagrange tensor product basis type \see OPENCMISS_BasisTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisExtendedLagrangeTPType = BASIS_EXTENDED_LAGRANGE_TP_TYPE !< Extendend Lagrange tensor product basis type \see OPENCMISS_BasisTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_BasisInterpolationSpecifications OPENCMISS::Basis::InterpolationSpecifications
  !> \brief Interpolation specification parameters
  !> \see OPENCMISS::BasisConstants,OPENCMISS
  !>@{ 
  INTEGER(INTG), PARAMETER :: CMISSBasisLinearLagrangeInterpolation = BASIS_LINEAR_LAGRANGE_INTERPOLATION !<Linear Lagrange interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisQuadraticLagrangeInterpolation = BASIS_QUADRATIC_LAGRANGE_INTERPOLATION !<Quadratic Lagrange interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisCubicLagrangeInterpolation = BASIS_CUBIC_LAGRANGE_INTERPOLATION !<Cubic Lagrange interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisCubicHermiteInterpolation = BASIS_CUBIC_HERMITE_INTERPOLATION !<Cubic Hermite interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisQuadratic1HermiteInterpolation = BASIS_QUADRATIC1_HERMITE_INTERPOLATION !<Quadratic Hermite (no derivative at xi=0) interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisQuadratic2HermiteInterpolation = BASIS_QUADRATIC2_HERMITE_INTERPOLATION !<Quadratic Hermite (no derivative at xi=1) interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisLinearSimplexInterpolation = BASIS_LINEAR_SIMPLEX_INTERPOLATION !<Linear Simplex interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisQuadraticSimplexInterpolation = BASIS_QUADRATIC_SIMPLEX_INTERPOLATION !<Quadratic Simplex interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisCubicSimplexInterpolation = BASIS_CUBIC_SIMPLEX_INTERPOLATION !<Cubic Simplex interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_BasisQuadratureTypes OPENCMISS::Basis::QuadratureTypes
  !> \brief Basis quadrature type parameters.
  !> \see OPENCMISS::BasisConstants,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSBasisGaussLegendreQuadrature = BASIS_GAUSS_LEGENDRE_QUADRATURE !<Gauss-Legendre quadrature \see OPENCMISS_BasisQuadratureTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisGaussLaguerreQuadrature = BASIS_GAUSS_LAGUERRE_QUADRATURE !<Gauss-Laguerre quadrature \see OPENCMISS_BasisQuadratureTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisGaussHermiteQuadrature = BASIS_GUASS_HERMITE_QUADRATURE !<Gauss-Hermite quadrature \see OPENCMISS_BasisQuadratureTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisAdaptiveGaussLegendreQuadrature = BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE !<Adaptive Gauss-Legendre quadrature \see OPENCMISS_BasisQuadratureTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisGaussSimplexQuadrature = BASIS_GAUSS_SIMPLEX_QUADRATURE !<Gauss-Legendre for Simplex elements quadrature \see OPENCMISS_BasisQuadratureTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_BasisXiCollapse OPENCMISS::Basis::XiCollapse
  !> \brief Basis Xi collapse parameters.
  !> \see OPENCMISS::Basis,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSBasisXiCollapsed = BASIS_XI_COLLAPSED !<The Xi direction is collapsed \see OPENCMISS_BasisXiCollapse,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisCollapsedAtXi0 = BASIS_COLLAPSED_AT_XI0 !<The Xi direction at the xi=0 end of this Xi direction is collapsed \see OPENCMISS_BasisXiCollapse,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisCollapsedAtXi1 = BASIS_COLLAPSED_AT_XI1 !<The Xi direction at the xi=1 end of this Xi direction is collapsed \see OPENCMISS_BasisXiCollapse,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBasisNotCollapsed = BASIS_NOT_COLLAPSED !<The Xi direction is not collapsed \see OPENCMISS_XiCollapse,OPENCMISS
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Returns the collapsed Xi flags for a basis.
  INTERFACE CMISSBasisCollapsedXiGet
    MODULE PROCEDURE CMISSBasisCollapsedXiGetNumber
    MODULE PROCEDURE CMISSBasisCollapsedXiGetObj
  END INTERFACE !CMISSBasisCollapsedXiGet
  
  !>Sets/changes the collapsed Xi flags for a basis.
  INTERFACE CMISSBasisCollapsedXiSet
    MODULE PROCEDURE CMISSBasisCollapsedXiSetNumber
    MODULE PROCEDURE CMISSBasisCollapsedXiSetObj
  END INTERFACE !CMISSBasisCollapsedXiSet
  
  !>Finishes the creation of a new basis. \see OPENCMISS::CMISSBasisCreateStart
  INTERFACE CMISSBasisCreateFinish
    MODULE PROCEDURE CMISSBasisCreateFinishNumber
    MODULE PROCEDURE CMISSBasisCreateFinishObj
  END INTERFACE !CMISSBasisCreateFinish
  
  !>Starts the creation of a new basis. \see OPENCMISS::CMISSBasisCreateFinish
  INTERFACE CMISSBasisCreateStart
    MODULE PROCEDURE CMISSBasisCreateFinishNumber
    MODULE PROCEDURE CMISSBasisCreateFinishObj
  END INTERFACE !CMISSBasisCreateFinish
  
  !>Destroys a basis.
  INTERFACE CMISSBasisDestroy
    MODULE PROCEDURE CMISSBasisDestroyNumber
    MODULE PROCEDURE CMISSBasisDestroyObj
  END INTERFACE !CMISSBasisDestroy

  !>Get the interpolation type in each Xi directions for a basis.
  INTERFACE CMISSBasisInterpolationXiGet
    MODULE PROCEDURE CMISSBasisInterpolationXiGetNumber
    MODULE PROCEDURE CMISSBasisInterpolationXiGetObj
  END INTERFACE !CMISSBasisInterpolationXiGet
  
  !>Sets/changes the interpolation type in each Xi directions for a basis.
  INTERFACE CMISSBasisInterpolationXiSet
    MODULE PROCEDURE CMISSBasisInterpolationXiSetNumber
    MODULE PROCEDURE CMISSBasisInterpolationXiSetObj
  END INTERFACE !CMISSBasisInterpolationXiSet
  
  !>Returns the number of local nodes in a basis.
  INTERFACE CMISSBasisNumberOfLocalNodesGet
    MODULE PROCEDURE CMISSBasisNumberOfLocalNodesGetNumber
    MODULE PROCEDURE CMISSBasisNumberOfLocalNodesGetObj
  END INTERFACE !CMISSBasisNumberOfLocalNodesGet
  
  !>Returns the number of Xi directions in a basis.
  INTERFACE CMISSBasisNumberOfXiGet
    MODULE PROCEDURE CMISSBasisNumberOfXiGetNumber
    MODULE PROCEDURE CMISSBasisNumberOfXiGetObj
  END INTERFACE !CMISSBasisNumberOfXiGet
  
  !>Sets/changes the number of Xi directions in a basis.
  INTERFACE CMISSBasisNumberOfXiSet
    MODULE PROCEDURE CMISSBasisNumberOfXiSetNumber
    MODULE PROCEDURE CMISSBasisNumberOfXiSetObj
  END INTERFACE !CMISSBasisNumberOfXiSet
  
  !>Returns the number of Gauss points in each Xi direction on a basis quadrature.
  INTERFACE CMISSBasisQuadratureNumberOfGaussXiGet
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiGetNumber
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiGetObj
  END INTERFACE !CMISSBasisQuadratureNumberOfGaussXiGet
  
  !>Sets/changes the number of Gauss points in each Xi direction on a basis quadrature.
  INTERFACE CMISSBasisQuadratureNumberOfGaussXiSet
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiSetNumber
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiSetObj
  END INTERFACE !CMISSBasisQuadratureNumberOfGaussXiSet
  
  !>Returns the order of quadrature for a basis quadrature.
  INTERFACE CMISSBasisQuadratureOrderGet
    MODULE PROCEDURE CMISSBasisQuadratureOrderGetNumber
    MODULE PROCEDURE CMISSBasisQuadratureOrderGetObj
  END INTERFACE !CMISSBasisQuadratureOrderGet
  
  !>Sets/changes the order of quadrature for a basis quadrature.
  INTERFACE CMISSBasisQuadratureOrderSet
    MODULE PROCEDURE CMISSBasisQuadratureOrderSetNumber
    MODULE PROCEDURE CMISSBasisQuadratureOrderSetObj
  END INTERFACE !CMISSBasisQuadratureOrderSet
  
  !>Returns the quadrature type for a basis quadrature.
  INTERFACE CMISSBasisQuadratureTypeGet
    MODULE PROCEDURE CMISSBasisQuadratureTypeGetNumber
    MODULE PROCEDURE CMISSBasisQuadratureTypeGetObj
  END INTERFACE !CMISSBasisQuadratureTypeGet
  
  !>Sets/changes the quadrature type for a basis quadrature.
  INTERFACE CMISSBasisQuadratureTypeSet
    MODULE PROCEDURE CMISSBasisQuadratureTypeSetNumber
    MODULE PROCEDURE CMISSBasisQuadratureTypeSetObj
  END INTERFACE !CMISSBasisQuadratureTypeSet
  
  !>Returns the type of a basis.
  INTERFACE CMISSBasisTypeGet
    MODULE PROCEDURE CMISSBasisTypeGetNumber
    MODULE PROCEDURE CMISSBasisTypeGetObj
  END INTERFACE !CMISSBasisTypeGet
  
  !>Sets/changes the type of a basis.
  INTERFACE CMISSBasisTypeSet
    MODULE PROCEDURE CMISSBasisTypeSetNumber
    MODULE PROCEDURE CMISSBasisTypeSetObj
  END INTERFACE !CMISSBasisTypeSet
  
  PUBLIC CMISSBasisLagrangeHermiteTPType,CMISSBasisSimplexType,CMISSBasisSerendipityType,CMISSBasisAuxilliaryType, &
    & CMISSBasisBSplineTPType,CMISSBasisFourierLagrangeHermiteTPType,CMISSBasisExtendedLagrangeTPType

  PUBLIC CMISSBasisLinearLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisCubicLagrangeInterpolation, &
    & CMISSBasisCubicHermiteInterpolation,CMISSBasisQuadratic1HermiteInterpolation,CMISSBasisQuadratic2HermiteInterpolation, &
    & CMISSBasisLinearSimplexInterpolation,CMISSBasisQuadraticSimplexInterpolation,CMISSBasisCubicSimplexInterpolation

  PUBLIC CMISSBasisGaussLegendreQuadrature,CMISSBasisGaussLaguerreQuadrature,CMISSBasisGaussHermiteQuadrature, &
    & CMISSBasisAdaptiveGaussLegendreQuadrature,CMISSBasisGaussSimplexQuadrature

  PUBLIC CMISSBasisXiCollapsed,CMISSBasisCollapsedAtXi0,CMISSBasisCollapsedAtXi1,CMISSBasisNotCollapsed

  PUBLIC CMISSBasisCollapsedXiGet,CMISSBasisCollapsedXiSet
  
  PUBLIC CMISSBasisCreateFinish,CMISSBasisCreateStart,CMISSBasisDestroy

  PUBLIC CMISSBasisInterpolationXiGet,CMISSBasisInterpolationXiSet

  PUBLIC CMISSBasisNumberOfLocalNodesGet

  PUBLIC CMISSBasisNumberOfXiGet,CMISSBasisNumberOfXiSet

  PUBLIC CMISSBasisQuadratureNumberOfGaussXiGet,CMISSBasisQuadratureNumberOfGaussXiSet

  PUBLIC CMISSBasisQuadratureOrderGet,CMISSBasisQuadratureOrderSet

  PUBLIC CMISSBasisQuadratureTypeGet,CMISSBasisQuadratureTypeSet
 
  PUBLIC CMISSBasisTypeGet,CMISSBasisTypeSet
 
!!==================================================================================================================================
!!
!! BOUNDARY_CONDITIONS_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters


  !> \addtogroup OPENCMISS_BoundaryConditionsConstants OPENCMISS::BoundaryConditions::Constants
  !> \brief Boundary conditions constants.
  !>@{  
  !> \addtogroup OPENCMISS_BoundaryConditionsTypes OPENCMISS::BoundaryConditions::Types
  !> \brief Boundary conditions type parameters.
  !> \see OPENCMISS::BoundaryConditions,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSBoundaryConditionNotFixed = BOUNDARY_CONDITION_NOT_FIXED !<The dof is not fixed. \see OPENCMISS_BoundaryConditionsTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBoundaryConditionFixed = BOUNDARY_CONDITION_FIXED !<The dof is fixed as a boundary condition. \see OPENCMISS_BoundaryConditionsTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSBoundaryConditionMixed = BOUNDARY_CONDITION_MIXED !<The dof is set as a mixed boundary condition. \see OPENCMISS_BoundaryConditionsTypes,OPENCMISS
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Destroys boundary conditions.
  INTERFACE CMISSBoundaryConditionsDestroy
    MODULE PROCEDURE CMISSBoundaryConditionsDestroyNumber
    MODULE PROCEDURE CMISSBoundaryConditionsDestroyObj
  END INTERFACE !CMISSBoundaryConditionsDestroy

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant.
  INTERFACE CMISSBoundaryConditionsAddConstant
    MODULE PROCEDURE CMISSBoundaryConditionsAddConstantNumber
    MODULE PROCEDURE CMISSBoundaryConditionsAddConstantObj
  END INTERFACE !CMISSBoundaryConditionsAddConstant

  !>Sets the value of the specified constant as a boundary condition on the specified constant.
  INTERFACE CMISSBoundaryConditionsSetConstant
    MODULE PROCEDURE CMISSBoundaryConditionsSetConstantNumber
    MODULE PROCEDURE CMISSBoundaryConditionsSetConstantObj
  END INTERFACE !CMISSBoundaryConditionsSetConstant

  !>Adds to the value of the element constant and sets this as a boundary condition on the specified element.
  INTERFACE CMISSBoundaryConditionsAddElement
    MODULE PROCEDURE CMISSBoundaryConditionsAddElementNumber
    MODULE PROCEDURE CMISSBoundaryConditionsAddElementObj
  END INTERFACE !CMISSBoundaryConditionsAddElement

  !>Sets the value of the specified element as a boundary condition on the specified element.
  INTERFACE CMISSBoundaryConditionsSetElement
    MODULE PROCEDURE CMISSBoundaryConditionsSetElementNumber
    MODULE PROCEDURE CMISSBoundaryConditionsSetElementObj
  END INTERFACE !CMISSBoundaryConditionsSetElement

  !>Adds to the value of the node constant and sets this as a boundary condition on the specified node.
  INTERFACE CMISSBoundaryConditionsAddNode
    MODULE PROCEDURE CMISSBoundaryConditionsAddNodeNumber
    MODULE PROCEDURE CMISSBoundaryConditionsAddNodeObj
  END INTERFACE !CMISSBoundaryConditionsAddNode

  !>Sets the value of the specified node as a boundary condition on the specified node.
  INTERFACE CMISSBoundaryConditionsSetNode
    MODULE PROCEDURE CMISSBoundaryConditionsSetNodeNumber
    MODULE PROCEDURE CMISSBoundaryConditionsSetNodeObj
  END INTERFACE !CMISSBoundaryConditionsSetNode

  !>Gets the boundary conditions for an equations set.
  INTERFACE CMISSEquationsSetBoundaryConditionsGet
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsGetNumber
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsGetObj
  END INTERFACE !CMISSEquationsSetBoundaryConditionsGet

  PUBLIC CMISSBoundaryConditionNotFixed,CMISSBoundaryConditionFixed,CMISSBoundaryConditionMixed

  PUBLIC CMISSBoundaryConditionsDestroy

  PUBLIC CMISSBoundaryConditionsAddConstant,CMISSBoundaryConditionsSetConstant

  PUBLIC CMISSBoundaryConditionsAddElement,CMISSBoundaryConditionsSetElement

  PUBLIC CMISSBoundaryConditionsAddNode,CMISSBoundaryConditionsSetNode

  PUBLIC CMISSEquationsSetBoundaryConditionsGet

!!==================================================================================================================================
!!
!! CONTROL_LOOP_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !> \addtogroup OPENCMISS_ControlLoopConstants OPENCMISS::ControlLoop::Constants
  !> \brief Control loops constants.
  !>@{  
  !> \addtogroup OPENCMISS_ControlLoopIdentifiers OPENCMISS::ControlLoop::Identifiers
  !> \brief The control loop identification parameters.
  !> \see OPENCMISS::ControlLoop,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSControlLoopNode = CONTROL_LOOP_NODE !<The identifier for a each "leaf" node in a control loop. \see OPENCMISS_ControlLoopIdentifiers,OPENCMISS
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Returns the current time parameters for a time control loop.
  INTERFACE CMISSControlLoopCurrentTimesGet
    MODULE PROCEDURE CMISSControlLoopCurrentTimesGetNumber0
    MODULE PROCEDURE CMISSControlLoopCurrentTimesGetNumber1
    MODULE PROCEDURE CMISSControlLoopCurrentTimesGetObj
  END INTERFACE !CMISSControlLoopCurrentTimesGet

  !>Destroy a control loop.
  INTERFACE CMISSControlLoopDestroy
    MODULE PROCEDURE CMISSControlLoopDestroyNumber0
    MODULE PROCEDURE CMISSControlLoopDestroyNumber1
    MODULE PROCEDURE CMISSControlLoopDestroyObj
  END INTERFACE !CMISSControlLoopDestroy

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  INTERFACE CMISSControlLoopGet
    MODULE PROCEDURE CMISSControlLoopGetNumber00
    MODULE PROCEDURE CMISSControlLoopGetNumber10
    MODULE PROCEDURE CMISSControlLoopGetNumber01
    MODULE PROCEDURE CMISSControlLoopGetNumber11
    MODULE PROCEDURE CMISSControlLoopGetObj0
    MODULE PROCEDURE CMISSControlLoopGetObj1
  END INTERFACE !CMISSControlLoopGet

  !>Sets/changes the iteration parameters for a fixed control loop. \todo need a get metod
  INTERFACE CMISSControlLoopIterationsSet
    MODULE PROCEDURE CMISSControlLoopIterationsSetNumber0
    MODULE PROCEDURE CMISSControlLoopIterationsSetNumber1
    MODULE PROCEDURE CMISSControlLoopIterationsSetObj
  END INTERFACE !CMISSControlLoopIterationsSet

  !>Sets/changes the maximum iterations for a while control loop. \todo need a get method
  INTERFACE CMISSControlLoopMaximumIterationsSet
    MODULE PROCEDURE CMISSControlLoopMaximumIterationsSetNumber0
    MODULE PROCEDURE CMISSControlLoopMaximumIterationsSetNumber1
    MODULE PROCEDURE CMISSControlLoopMaximumIterationsSetObj
  END INTERFACE !CMISSControlLoopMaximumIterationsSet

  !>Returns the number of sub loops for a control loop.
  INTERFACE CMISSControlLoopNumberOfSubLoopsGet
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsGetNumber0
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsGetNumber1
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsGetObj
  END INTERFACE !CMISSControlLoopNumberOfSubLoopsGet

  !>Sets/changes the number of sub loops for a control loop. \todo is this really a public method???
  INTERFACE CMISSControlLoopNumberOfSubLoopsSet
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsSetNumber0
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsSetNumber1
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsSetObj
  END INTERFACE !CMISSControlLoopNumberOfSubLoopsGet

  !>Returns the time parameters for a time control loop.
  INTERFACE CMISSControlLoopTimesGet
    MODULE PROCEDURE CMISSControlLoopTimesGetNumber0
    MODULE PROCEDURE CMISSControlLoopTimesGetNumber1
    MODULE PROCEDURE CMISSControlLoopTimesGetObj
  END INTERFACE !CMISSControlLoopTimesGet

  !>Sets/Changes the time parameters for a time control loop.
  INTERFACE CMISSControlLoopTimesSet
    MODULE PROCEDURE CMISSControlLoopTimesSetNumber0
    MODULE PROCEDURE CMISSControlLoopTimesSetNumber1
    MODULE PROCEDURE CMISSControlLoopTimesSetObj
  END INTERFACE !CMISSControlLoopTimesSet

  !>Sets/Changes the loop type for a control loop. \todo Is this really a public method? \todo need a get method
  INTERFACE CMISSControlLoopTypeSet
    MODULE PROCEDURE CMISSControlLoopTypeSetNumber0
    MODULE PROCEDURE CMISSControlLoopTypeSetNumber1
    MODULE PROCEDURE CMISSControlLoopTypeSetObj
  END INTERFACE !CMISSControlLoopTypeSet

   PUBLIC CMISSControlLoopCurrentTimesGet
   
  PUBLIC CMISSControlLoopDestroy

  PUBLIC CMISSControlLoopGet

  PUBLIC CMISSControlLoopIterationsSet

  PUBLIC CMISSControlLoopMaximumIterationsSet

  PUBLIC CMISSControlLoopNumberOfSubLoopsGet,CMISSControlLoopNumberOfSubLoopsSet

  PUBLIC CMISSControlLoopTimesGet,CMISSControlLoopTimesSet

  PUBLIC CMISSControlLoopTypeSet
 

!!==================================================================================================================================
!!
!! COORDINATE_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !> \addtogroup OPENCMISS_CoordinateConstants OPENCMISS::Coordinate::Constants
  !> \brief Coordinate constants.
  !>@{  
  !> \addtogroup OPENCMISS_CoordinateSystemTypes OPENCMISS::Coordinate::SystemTypes
  !> \brief Coordinate system type parameters.
  !> \see OPENCMISS::Coordinate,OPENCMISS
  !>@{ 
  INTEGER(INTG), PARAMETER :: CMISSCoordinateRectangularCartesianType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE !<Rectangular Cartesian coordinate system type \see OPENCMISS_CoordinateSystemTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateCylindricalPolarType = COORDINATE_CYLINDRICAL_POLAR_TYPE !<Cylindrical polar coordinate system type \see OPENCMISS_CoordinateSystemTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateSphericalPolarType = COORDINATE_SPHERICAL_POLAR_TYPE !<Spherical polar coordinate system type \see OPENCMISS_CoordinateSystemTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateProlateSpheroidalType = COORDINATE_PROLATE_SPHEROIDAL_TYPE !<Prolate spheroidal coordinate system type \see OPENCMISS_CoordinateSystemTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateOblateSpheroidalType = COORDINATE_OBLATE_SPHEROIDAL_TYPE !<Oblate spheroidal coordinate system type \see OPENCMISS_CoordinateSystemTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_CoordinateRadialInterpolations OPENCMISS::Coordinate::RadialInterpolations
  !> \brief The type of radial interpolation for polar coordinate systems
  !> \see OPENCMISS::Coordinate,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSCoordinateNoRadialInterpolationType = COORDINATE_NO_RADIAL_INTERPOLATION_TYPE !<No radial interpolation \see OPENCMISS_CoordinateRadialInterpolations,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateRadialInterpolationType = COORDINATE_RADIAL_INTERPOLATION_TYPE !<r radial interpolation \see OPENCMISS_CoordinateRadialInterpolations,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateRadialSquaredInterpolationType = COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE !<r^2 radial interpolation \see OPENCMISS_CoordinateRadialInterpolations,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSCoordinateRadialCubedInterpolationType = COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE !<r^3 radial interpolation \see OPENCMISS_CoordinateRadialInterpolations,OPENCMISS
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Finishes the creation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemCreateStart
  INTERFACE CMISSCoordinateSystemCreateFinish
    MODULE PROCEDURE CMISSCoordinateSystemCreateFinishNumber
    MODULE PROCEDURE CMISSCoordinateSystemCreateFinishObj
  END INTERFACE !CMISSCoordinateSystemCreateFinish

  !>Starts the creation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemCreateFinish
  INTERFACE CMISSCoordinateSystemCreateStart
    MODULE PROCEDURE CMISSCoordinateSystemCreateStartNumber
    MODULE PROCEDURE CMISSCoordinateSystemCreateStartObj
  END INTERFACE !CMISSCoordinateSystemCreateStart

  !>Destorys a coordinate system.
  INTERFACE CMISSCoordinateSystemDestroy
    MODULE PROCEDURE CMISSCoordinateSystemDestroyNumber
    MODULE PROCEDURE CMISSCoordinateSystemDestroyObj
  END INTERFACE !CMISSCoordinateSystemDestroy

  !>Returns the coordinate system dimension. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemDimensionGet
    MODULE PROCEDURE CMISSCoordinateSystemDimensionGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemDimensionGetObj
  END INTERFACE !CMISSCoordinateSystemDimensionGet

  !>Sets/changes the coordinate system dimension. \todo fix pointers
  INTERFACE CMISSCoordinateSystemDimensionSet
    MODULE PROCEDURE CMISSCoordinateSystemDimensionSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemDimensionSetObj
  END INTERFACE !CMISSCoordinateSystemDimensionSet

  !>Returns the coordinate system focus. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemFocusGet
    MODULE PROCEDURE CMISSCoordinateSystemFocusGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemFocusGetObj
  END INTERFACE !CMISSCoordinateSystemFocusGet
    
  !>Sets/changes the coordinate system focus. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemFocusSet
    MODULE PROCEDURE CMISSCoordinateSystemFocusSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemFocusSetObj
  END INTERFACE !CMISSCoordinateSystemFocusSet

  !>Returns the coordinate system radial interpolation type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemRadialInterpolationGet
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationGetObj
  END INTERFACE !CMISSCoordinateSystemRadialInterpolationGet
    
  !>Sets/changes the coordinate system radial interpolation type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemRadialInterpolationSet
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationSetObj
  END INTERFACE !CMISSCoordinateSystemRadialInterpolationSet
    
  !>Returns the coordinate system type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemTypeGet
    MODULE PROCEDURE CMISSCoordinateSystemTypeGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemTypeGetObj
  END INTERFACE !CMISSCoordinateSystemTypeGet
    
  !>Sets/changes the coordinate system type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemTypeSet
    MODULE PROCEDURE CMISSCoordinateSystemTypeSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemTypeSetObj
  END INTERFACE !CMISSCoordinateSystemTypeSet

  !>Returns the coordinate system orign. 
  INTERFACE CMISSCoordinateSystemOriginGet
    MODULE PROCEDURE CMISSCoordinateSystemOriginGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOriginGetObj
  END INTERFACE !CMISSCoordinateSystemOriginGet

  !>Sets/changes the coordinate system orign. 
  INTERFACE CMISSCoordinateSystemOriginSet
    MODULE PROCEDURE CMISSCoordinateSystemOriginSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOriginSetObj
  END INTERFACE !CMISSCoordinateSystemOriginSet

  !>Returns the coordinate system orientation. 
  INTERFACE CMISSCoordinateSystemOrientationGet
    MODULE PROCEDURE CMISSCoordinateSystemOrientationGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOrientationGetObj
  END INTERFACE !CMISSCoordinateSystemOrientationGet

  !>Sets/changes the coordinate system orientation. 
  INTERFACE CMISSCoordinateSystemOrientationSet
    MODULE PROCEDURE CMISSCoordinateSystemOrientationSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOrientationSetObj
  END INTERFACE !CMISSCoordinateSystemOrientationSet

  PUBLIC CMISSCoordinateRectangularCartesianType,CMISSCoordinateCylindricalPolarType,CMISSCoordinateSphericalPolarType, &
    & CMISSCoordinateProlateSpheroidalType,CMISSCoordinateOblateSpheroidalType

  PUBLIC CMISSCoordinateNoRadialInterpolationType,CMISSCoordinateRadialInterpolationType, &
    & CMISSCoordinateRadialSquaredInterpolationType,CMISSCoordinateRadialCubedInterpolationType

  PUBLIC CMISSCoordinateSystemCreateFinish,CMISSCoordinateSystemCreateStart

  PUBLIC CMISSCoordinateSystemDestroy
  
  PUBLIC CMISSCoordinateSystemDimensionGet,CMISSCoordinateSystemDimensionSet

  PUBLIC CMISSCoordinateSystemFocusGet,CMISSCoordinateSystemFocusSet

  PUBLIC CMISSCoordinateSystemRadialInterpolationGet,CMISSCoordinateSystemRadialInterpolationSet

  PUBLIC CMISSCoordinateSystemTypeGet,CMISSCoordinateSystemTypeSet

  PUBLIC CMISSCoordinateSystemOriginGet,CMISSCoordinateSystemOriginSet

  PUBLIC CMISSCoordinateSystemOrientationGet,CMISSCoordinateSystemOrientationSet
  
!!==================================================================================================================================
!!
!! EQUATIONS_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters
  
  !> \addtogroup OPENCMISS_EquationsConstants OPENCMISS::Equations::Constants
  !> \brief Equations  constants.
  !>@{
  !> \addtogroup OPENCMISS_EquationsOutputTypes OPENCMISS::Equations::OutputTypes
  !> \brief Equations output types
  !> \see OPENCMISS::Equations,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsNoOutput = EQUATIONS_NO_OUTPUT!<No output from the equations \see OPENCMISS_EquationsOutputTypes,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsTimingOutput = EQUATIONS_TIMING_OUTPUT !<Timing information output. \see OPENCMISS_EquationsOutputTypes,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsMatrixOutput = EQUATIONS_MATRIX_OUTPUT !<All below and equation matrices output. \see OPENCMISS_EquationsOutputTypes,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsElementMatrixOutput = EQUATIONS_ELEMENT_MATRIX_OUTPUT !<All below and element matrices output. \see OPENCMISS_EquationsOutputTypes,OPENCMISS   
  !>@}
  !> \addtogroup OPENCMISS_EquationsSparsityTypes OPENCMISS::Equations::SparsityTypes
  !> \brief Equations sparsity types
  !> \see OPENCMISS::Equations,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsSparseMatrices = EQUATIONS_SPARSE_MATRICES !<Use sparse matrices for the equations. \see OPENCMISS_EquationsSparsityTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsFullMatrices = EQUATIONS_FULL_MATRICES !<Use fully populated matrices for the equations. \see OPENCMISS_EquationsSparsityTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_EquationsLumpingTypes OPENCMISS::Equations::LumpingTypes
  !> \brief Equations  lumping types
  !> \see OPENCMISS_EquationsSparsityTypes,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsUnlumpedMatrices = EQUATIONS_UNLUMPED_MATRICES !<The equations matrices are not lumped. \see OPENCMISS_EquationsLumpingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsLumpedMatrices = EQUATIONS_LUMPED_MATRICES !<The equations matrices are "mass" lumped. \see OPENCMISS_EquationsLumpingTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_EquationsLinearityTypes OPENCMISS::Equations::LinearityTypes
  !> \brief The equations linearity types
  !> \see OPENCMISS::Equations,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsLinear = EQUATIONS_LINEAR !<The equations are linear. \see OPENCMISS_EquationsLinearityTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsNonlinear = EQUATIONS_NONLINEAR !<The equations are non-linear. \see \see OPENCMISS_EquationsLinearityTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsNonlinearBCs = EQUATIONS_NONLINEAR_BCS !<The equations have non-linear boundary conditions. \see \see OPENCMISS_EquationsLinearityTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_EquationsTimeDepedenceTypes OPENCMISS::Equations::TimeDepedenceTypes
  !> \brief The equations time dependence types
  !> \see OPENCMISS::Equations,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsStatic = EQUATIONS_STATIC !<The equations are static and have no time dependence. \see OPENCMISS_EquationsTimeDepedenceTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsQuasistatic = EQUATIONS_QUASISTATIC !<The equations are quasi-static. \see OPENCMISS_EquationsTimeDepedenceTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsFirstOrderDynamic = EQUATIONS_FIRST_ORDER_DYNAMIC !<The equations are first order dynamic. \see OPENCMISS_EquationsTimeDepedenceTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSecondOrderDynamic = EQUATIONS_SECOND_ORDER_DYNAMIC !<The equations are a second order dynamic. \see OPENCMISS_EquationsTimeDepedenceTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsTimeStepping = EQUATIONS_TIME_STEPPING !<The equations are for time stepping. \see OPENCMISS_EquationsTimeDepedenceTypes,OPENCMISS
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Destroys equations for an equations set.
  INTERFACE CMISSEquationsDestroy
    MODULE PROCEDURE CMISSEquationsDestroyNumber
    MODULE PROCEDURE CMISSEquationsDestroyObj
  END INTERFACE !CMISSEquationsDestroy

  !>Gets the linearity type for equations.
  INTERFACE CMISSEquationsLinearityTypeGet
    MODULE PROCEDURE CMISSEquationsLinearityTypeGetNumber
    MODULE PROCEDURE CMISSEquationsLinearityTypeGetObj
  END INTERFACE !CMISSEquationsLinearityTypeGet

  !>Gets the lumping type for equations.
  INTERFACE CMISSEquationsLumpingTypeGet
    MODULE PROCEDURE CMISSEquationsLumpingTypeGetNumber
    MODULE PROCEDURE CMISSEquationsLumpingTypeGetObj
  END INTERFACE !CMISSEquationsLumpingTypeGet

  !>Sets/changes the lumping type for equations.
  INTERFACE CMISSEquationsLumpingTypeSet
    MODULE PROCEDURE CMISSEquationsLumpingTypeSetNumber
    MODULE PROCEDURE CMISSEquationsLumpingTypeSetObj
  END INTERFACE !CMISSEquationsLumpingTypeSet

  !>Gets the output type for equations.
  INTERFACE CMISSEquationsOutputTypeGet
    MODULE PROCEDURE CMISSEquationsOutputTypeGetNumber
    MODULE PROCEDURE CMISSEquationsOutputTypeGetObj
  END INTERFACE !CMISSEquationsOutputTypeGet

  !>Sets/changes the output type for equations.
  INTERFACE CMISSEquationsOutputTypeSet
    MODULE PROCEDURE CMISSEquationsOutputTypeSetNumber
    MODULE PROCEDURE CMISSEquationsOutputTypeSetObj
  END INTERFACE !CMISSEquationsOutputTypeSet

 !>Gets the sparsity type for equations.
  INTERFACE CMISSEquationsSparsityTypeGet
    MODULE PROCEDURE CMISSEquationsSparsityTypeGetNumber
    MODULE PROCEDURE CMISSEquationsSparsityTypeGetObj
  END INTERFACE !CMISSEquationsSparsityTypeGet

  !>Sets/changes the sparsity type for equations.
  INTERFACE CMISSEquationsSparsityTypeSet
    MODULE PROCEDURE CMISSEquationsSparsityTypeSetNumber
    MODULE PROCEDURE CMISSEquationsSparsityTypeSetObj
  END INTERFACE !CMISSEquationsSparsityTypeSet

  !>Gets the time dependence type for equations.
  INTERFACE CMISSEquationsTimeDependenceTypeGet
    MODULE PROCEDURE CMISSEquationsTimeDependenceTypeGetNumber
    MODULE PROCEDURE CMISSEquationsTimeDependenceTypeGetObj
  END INTERFACE !CMISSEquationsTimeDependenceTypeGet

 PUBLIC CMISSEquationsNoOutput,CMISSEquationsTimingOutput,CMISSEquationsMatrixOutput,CMISSEquationsElementMatrixOutput

  PUBLIC CMISSEquationsSparseMatrices,CMISSEquationsFullMatrices

  PUBLIC CMISSEquationsUnlumpedMatrices,CMISSEquationsLumpedMatrices

  PUBLIC CMISSEquationsLinear,CMISSEquationsNonlinear,CMISSEquationsNonlinearBCs

  PUBLIC CMISSEquationsStatic,CMISSEquationsQuasistatic,CMISSEquationsFirstOrderDynamic,CMISSEquationsSecondOrderDynamic, &
    & CMISSEquationsTimeStepping

  PUBLIC CMISSEquationsDestroy

  PUBLIC CMISSEquationsLinearityTypeGet

  PUBLIC CMISSEquationsLumpingTypeGet,CMISSEquationsLumpingTypeSet

  PUBLIC CMISSEquationsOutputTypeGet,CMISSEquationsOutputTypeSet

  PUBLIC CMISSEquationsSparsityTypeGet,CMISSEquationsSparsityTypeSet

  PUBLIC CMISSEquationsTimeDependenceTypeGet

!!==================================================================================================================================
!!
!! EQUATIONS_SET_CONSTANTS
!!
!!==================================================================================================================================

  !Module parameters
  
  !> \addtogroup OPENCMISS_EquationsSetConstants OPENCMISS::EquationsSet::Constants
  !> \brief Equations set constants.
  !>@{
  !> \addtogroup OPENCMISS_EquationsSetClasses OPENCMISS::EquationsSet::Classes
  !> \brief Equations set classes.
  !> \see OPENCMISS::EquationsSet,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetNoClass = EQUATIONS_SET_NO_CLASS !<No equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetElasticityClass = EQUATIONS_SET_ELASTICITY_CLASS !<Elasticity equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFluidMechanicsClass = EQUATIONS_SET_FLUID_MECHANICS_CLASS !<Fluid Mechanics equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetElectroMechanicsClass = EQUATIONS_SET_ELECTROMAGNETICS_CLASS !<Electromagnetics equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetClassicalFieldClass = EQUATIONS_SET_CLASSICAL_FIELD_CLASS !<Classical Field equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS   
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetBioelectricsClass = EQUATIONS_SET_BIOELECTRICS_CLASS !<Bioelectrics equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS     
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetModalClass = EQUATIONS_SET_MODAL_CLASS !<Modal equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS     
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFittingClass = EQUATIONS_SET_FITTING_CLASS !<Fitting equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS     
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetOptimisationClass = EQUATIONS_SET_OPTIMISATION_CLASS !<Optimisation equations set class \see OPENCMISS_EquationsSetClasses,OPENCMISS     
  !>@}
  !> \addtogroup OPENCMISS_EquationsSetTypes OPENCMISS::EquationsSet::Types
  !> \brief Equations set Types.
  !> \see OPENCMISS::EquationsSet,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetNoType = EQUATIONS_SET_NO_TYPE !<No equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLinearElasticityType = EQUATIONS_SET_LINEAR_ELASTICITY_TYPE !<Linear elasticity equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFiniteElasticityType = EQUATIONS_SET_FINITE_ELASTICITY_TYPE !<Finite elasticity equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetStokesEquationType = EQUATIONS_SET_STOKES_EQUATION_TYPE !<Stokes equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetNavierStokesType = EQUATIONS_SET_NAVIER_STOKES_TYPE !<Navier-Stokes equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetDarcyEquationType = EQUATIONS_SET_DARCY_EQUATION_TYPE !<Darcy equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetElectrostaticType = EQUATIONS_SET_ELECTROSTATIC_TYPE !<Electrostatic equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetMagnetoStaticType = EQUATIONS_SET_MAGNETOSTATIC_TYPE !<Magnetostatic equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetMaxwellsEquationType = EQUATIONS_SET_MAXWELLS_EQUATIONS_TYPE !<Maxwells equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceEquationType = EQUATIONS_SET_LAPLACE_EQUATION_TYPE !<Laplace equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetPoissonEquationType = EQUATIONS_SET_POISSON_EQUATION_TYPE !<Poisson equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetHelmholtzEquationType = EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE !<Helmholtz equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetWaveEquationType = EQUATIONS_SET_WAVE_EQUATION_TYPE !<Wave equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetDiffusionEquationType = EQUATIONS_SET_DIFFUSION_EQUATION_TYPE !<Diffusion equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetAdvectionDiffusionEquationType = EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE !<Advection-Diffusion equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetReactionDiffusionEquationType = EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE !<Reaction-Diffusion equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetBiharmonicEquationType = EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE !<Biharmonic equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetMonodomainEquationType = EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE !<Monodomain equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetBidomainEquationType = EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE !<Bidomain equation equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLinearElasticModalType = EQUATIONS_SET_LINEAR_ELASTIC_MODAL_TYPE !<Linear elasticity modal equations set type \see OPENCMISS_EquationsSetTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_EquationsSetSubtypes OPENCMISS::EquationsSet::Subtypes
  !> \brief Equations set subtypes.
  !> \see OPENCMISS::EquationsSet,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetNoSubtype = EQUATIONS_SET_NO_SUBTYPE !<No equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetThreeDimensionalLinearElasticitySubtype = EQUATIONS_SET_THREE_DIMENSIONAL_LINEAR_ELASTICITY_SUBTYPE !<Three dimensional linear elasticity equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetPlaneStressSubtype = EQUATIONS_SET_PLANE_STRESS_SUBTYPE !<Plane stress linear elasticity equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetPlaneStrainSubtype = EQUATIONS_SET_PLANE_STRAIN_SUBTYPE !<Plane strain linear elasticity equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetPlateSubtype = EQUATIONS_SET_PLATE_SUBTYPE !<Plate linear elasticity equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetShellSubtype = EQUATIONS_SET_SHELL_SUBTYPE !<Shell linear elasticity equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetStaticStokesSubtype = EQUATIONS_SET_STATIC_STOKES_SUBTYPE !<Static Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceStokesSubtype = EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE !<Laplace type Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetTransientStokesSubtype = EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE !<Transient Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetOptimisedStokesSubtype = EQUATIONS_SET_OPTIMISED_STOKES_SUBTYPE !<Optimised Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetStaticNavierStokesSubtype = EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE !<Static Navier-Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceNavierStokesSubtype = EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE !<Laplace type Navier-Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetTransientNavierStokesSubtype = EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE !<Transient Navier-Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetOptimisedNavierStokesSubtype = EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE !<Optimised Navier-Stokes equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetStandardDarcySubtype = EQUATIONS_SET_STANDARD_DARCY_SUBTYPE !<Standard Darcy equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetStandardLaplaceSubtype = EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE !<Standard Laplace equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetGeneralisedLaplaceSubtype = EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE !<Generalised Laplace equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetConstantSourcePoissonSubtype = EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE !<Constant source Poisson equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLinearSourcePoissonSubtype = EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE !<Linear source Poisson equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetQuadraticSourcePoissonSubtype = EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE !<Quadratic source Poisson equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetExponentialSourcePoissonSubtype = EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE !<Exponential source Poisson equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetNoSourceHelmholtzSubtype = EQUATIONS_SET_NO_SOURCE_HELMHOLTZ_SUBTYPE !<No source Helmholtz equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetNoSourceDiffusionSubtype = EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE !<No source diffusion equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFirstBidomainSubtype = EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE !<First bidomain equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetSecondBidomainSubtype = EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE !<Second bidomain equations set subtype \see OPENCMISS_EquationsSetSubtypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_EquationsSetSolutionMethods OPENCMISS::EquationsSet::SolutionMethods
  !> \brief The solution method parameters
  !> \see OPENCMISS::EquationsSet,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFEMSolutionMethod = EQUATIONS_SET_FEM_SOLUTION_METHOD !<Finite Element Method solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetBEMSolutionMethod = EQUATIONS_SET_BEM_SOLUTION_METHOD !<Boundary Element Method solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFDSolutionMethod = EQUATIONS_SET_FD_SOLUTION_METHOD !<Finite Difference solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetFVSolutionMethod = EQUATIONS_SET_FV_SOLUTION_METHOD !<Finite Volume solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetGFEMSolutionMethod = EQUATIONS_SET_GFEM_SOLUTION_METHOD !<Grid-based Finite Element Method solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetGFDSolutionMethod = EQUATIONS_SET_GFD_SOLUTION_METHOD !<Grid-based Finite Difference solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetGFVSolutionMethod = EQUATIONS_SET_GFV_SOLUTION_METHOD !<Grid-based Finite Volume solution method. \see OPENCMISS_EquationsSetSolutionMethods,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_EquationsSetAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes
  !> \brief The analytic function types.
  !> \see OPENCMISS::EquationsSet,OPENCMISS
  !>@{
  !> \addtogroup OPENCMISS_EquationsSetLaplaceAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::Laplace
  !> \brief The analytic function types for a Laplace equation
  !> \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
  !>@{  
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceEquationTwoDim1 = EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1 !<u=x**2+2*x*y-y**2 \see OPENCMISS_EquationsSetLaplaceAnalyticFunctionTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceEquationTwoDim2 = EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2 !<u=cos(x)cosh(y) \see OPENCMISS_EquationsSetLaplaceAnalyticFunctionTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceEquationThreeDim1 = EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1 !<u=x**2-2*y**2+z**2 \see OPENCMISS_EquationsSetLaplaceAnalyticFunctionTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetLaplaceEquationThreeDim2 = EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2 !<u=cos(x)*cosh(y)*z \see OPENCMISS_EquationsSetLaplaceAnalyticFunctionTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_PoissonAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::Poisson
  !> \brief The analytic function types for a Poisson equation.
  !> \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
  !>@{  
  INTEGER(INTG), PARAMETER :: CMISSEquationsSetPoissonTwoDim1 = EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1 !<u=ln(4/(x+y+1^2)) \see OPENCMISS_EquationsSetPoissonAnalyticFunctionTypes,OPENCMISS
  !>@}
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSEquationsSetNoClass,CMISSEquationsSetElasticityClass,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetElectroMechanicsClass,CMISSEquationsSetClassicalFieldClass,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetModalClass,CMISSEquationsSetFittingClass,CMISSEquationsSetOptimisationClass

  PUBLIC CMISSEquationsSetNoType,CMISSEquationsSetLinearElasticityType,CMISSEquationsSetFiniteElasticityType, &
    & CMISSEquationsSetStokesEquationType,CMISSEquationsSetNavierStokesType,CMISSEquationsSetDarcyEquationType, &
    & CMISSEquationsSetElectrostaticType,CMISSEquationsSetMagnetoStaticType,CMISSEquationsSetMaxwellsEquationType, &
    & CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetPoissonEquationType,CMISSEquationsSetHelmholtzEquationType, &
    & CMISSEquationsSetWaveEquationType,CMISSEquationsSetDiffusionEquationType,CMISSEquationsSetAdvectionDiffusionEquationType, &
    & CMISSEquationsSetReactionDiffusionEquationType,CMISSEquationsSetBiharmonicEquationType, &
    & CMISSEquationsSetMonodomainEquationType,CMISSEquationsSetBidomainEquationType,CMISSEquationsSetLinearElasticModalType

  PUBLIC CMISSEquationsSetNoSubtype,CMISSEquationsSetThreeDimensionalLinearElasticitySubtype,CMISSEquationsSetPlaneStressSubtype, &
    & CMISSEquationsSetPlaneStrainSubtype,CMISSEquationsSetPlateSubtype,CMISSEquationsSetShellSubtype, &
    & CMISSEquationsSetStaticStokesSubtype,CMISSEquationsSetLaplaceStokesSubtype,CMISSEquationsSetTransientStokesSubtype, &
    & CMISSEquationsSetOptimisedStokesSubtype,CMISSEquationsSetStaticNavierStokesSubtype, &
    & CMISSEquationsSetLaplaceNavierStokesSubtype,CMISSEquationsSetTransientNavierStokesSubtype,&
    & CMISSEquationsSetOptimisedNavierStokesSubtype,CMISSEquationsSetStandardDarcySubtype,CMISSEquationsSetStandardLaplaceSubtype, &
    & CMISSEquationsSetGeneralisedLaplaceSubtype,CMISSEquationsSetConstantSourcePoissonSubtype, &
    & CMISSEquationsSetLinearSourcePoissonSubtype,CMISSEquationsSetQuadraticSourcePoissonSubtype, &
    & CMISSEquationsSetExponentialSourcePoissonSubtype,CMISSEquationsSetNoSourceHelmholtzSubtype, &
    & CMISSEquationsSetNoSourceDiffusionSubtype,CMISSEquationsSetFirstBidomainSubtype,CMISSEquationsSetSecondBidomainSubtype

  PUBLIC CMISSEquationsSetFEMSolutionMethod,CMISSEquationsSetBEMSolutionMethod,CMISSEquationsSetFDSolutionMethod, &
    & CMISSEquationsSetFVSolutionMethod,CMISSEquationsSetGFEMSolutionMethod,CMISSEquationsSetGFDSolutionMethod, &
    & CMISSEquationsSetGFVSolutionMethod

  PUBLIC CMISSEquationsSetLaplaceEquationTwoDim1,CMISSEquationsSetLaplaceEquationTwoDim2, &
    & CMISSEquationsSetLaplaceEquationThreeDim1,CMISSEquationsSetLaplaceEquationThreeDim2

  PUBLIC CMISSEquationsSetPoissonTwoDim1
  
!!==================================================================================================================================
!!
!! EQUATIONS_SET_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters
  
  !Module types

  !Module variables

  !Interfaces

  !>Finish the creation of a analytic solution for an equations set. \see OPENCMISS::CMISSEquationsSetAnalyticCreateStart
  INTERFACE CMISSEquationsSetAnalyticCreateFinish
    MODULE PROCEDURE CMISSEquationsSetAnalyticCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetAnalyticCreateFinishObj
  END INTERFACE !CMISSEquationsSetAnalyticCreateFinish
  
  !>Start the creation of a analytic solution for an equations set. \see OPENCMISS::CMISSEquationsSetAnalyticCreateFinish
  INTERFACE CMISSEquationsSetAnalyticCreateStart
    MODULE PROCEDURE CMISSEquationsSetAnalyticCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetAnalyticCreateStartObj
  END INTERFACE !CMISSEquationsSetAnalyticCreateStart
  
  !>Destroy the analytic solution for an equations set.
  INTERFACE CMISSEquationsSetAnalyticDestroy
    MODULE PROCEDURE CMISSEquationsSetAnalyticDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetAnalyticDestroyObj
  END INTERFACE !CMISSEquationsSetAnalyticDestroy
  
  !>Set boundary conditions for an equation set according to the analytic equations.
  INTERFACE CMISSEquationsSetBoundaryConditionsAnalytic
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsAnalyticNumber
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsAnalyticObj
  END INTERFACE !CMISSEquationsSetBoundaryConditionsAnalytic
  
  !>Finish the creation of boundary conditions for an equation set. \see OPENCMISS::CMISSEquationsSetBoundaryConditionsCreateStart
  INTERFACE CMISSEquationsSetBoundaryConditionsCreateFinish
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsCreateFinishObj
  END INTERFACE !CMISSEquationsSetBoundaryConditionsCreateFinish
  
  !>Start the creation of boundary conditions for an equation set. \see OPENCMISS::CMISSEquationsSetBoundaryConditionsCreateFinish
  INTERFACE CMISSEquationsSetBoundaryConditionsCreateStart
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsCreateStartObj
  END INTERFACE !CMISSEquationsSetBoundaryConditionsCreateStart
  
  !>Destroy the boundary conditions for an equations set.
  INTERFACE CMISSEquationsSetBoundaryConditionsDestroy
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsDestroyObj
  END INTERFACE !CMISSEquationsSetBoundaryConditionsDestroy
  
  !>Finish the creation of an equations set. \see OPENCMISS::CMISSEquationsSetCreateStart
  INTERFACE CMISSEquationsSetCreateFinish
    MODULE PROCEDURE CMISSEquationsSetCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetCreateFinishObj
  END INTERFACE !CMISSEquationsSetCreateFinish
  
  !>Start the creation of an equations set on a region. \see OPENCMISS::CMISSEquationsSetCreateFinish
  INTERFACE CMISSEquationsSetCreateStart
    MODULE PROCEDURE CMISSEquationsSetCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetCreateStartObj
  END INTERFACE !CMISSEquationsSetCreateStart
  
  !>Destroy an equations set. 
  INTERFACE CMISSEquationsSetDestroy
    MODULE PROCEDURE CMISSEquationsSetDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetDestroyObj
  END INTERFACE !CMISSEquationsSetDestroy
  
  !>Finish the creation of dependent variables for an equations set. \see OPENCMISS::CMISSEquationsSetDependentCreateStart
  INTERFACE CMISSEquationsSetDependentCreateFinish
    MODULE PROCEDURE CMISSEquationsSetDependentCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetDependentCreateFinishObj
  END INTERFACE !CMISSEquationsSetDependentCreateFinish
  
  !>Start the creation of dependent variables for an equations set. \see OPENCMISS::CMISSEquationsSetDependentCreateFinish
  INTERFACE CMISSEquationsSetDependentCreateStart
    MODULE PROCEDURE CMISSEquationsSetDependentCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetDependentCreateStartObj
  END INTERFACE !CMISSEquationsSetDependentCreateStart
  
  !>Destroy the dependent variables for an equations set.
  INTERFACE CMISSEquationsSetDependentDestroy
    MODULE PROCEDURE CMISSEquationsSetDependentDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetDependentDestroyObj
  END INTERFACE !CMISSEquationsSetDependentDestroy
  
  !>Finish the creation of equations for an equations set. \see OPENCMISS::CMISSEquationsSetEquationsCreateStart
  INTERFACE CMISSEquationsSetEquationsCreateFinish
    MODULE PROCEDURE CMISSEquationsSetEquationsCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetEquationsCreateFinishObj
  END INTERFACE !CMISSEquationsSetEquationsCreateFinish
  
  !>Start the creation of equations for an equations set. \see OPENCMISS::CMISSEquationsSetEquationsCreateFinish
  INTERFACE CMISSEquationsSetEquationsCreateStart
    MODULE PROCEDURE CMISSEquationsSetEquationsCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetEquationsCreateStartObj
  END INTERFACE !CMISSEquationsSetEquationsCreateStart
  
  !>Destroy the equations for an equations set.
  INTERFACE CMISSEquationsSetEquationsDestroy
    MODULE PROCEDURE CMISSEquationsSetEquationsDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetEquationsDestroyObj
  END INTERFACE !CMISSEquationsSetEquationsDestroy
  
  !>Finish the creation of materials for an equations set. \see OPENCMISS::CMISSEquationsSetMaterialsCreateStart
  INTERFACE CMISSEquationsSetMaterialsCreateFinish
    MODULE PROCEDURE CMISSEquationsSetMaterialsCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetMaterialsCreateFinishObj
  END INTERFACE !CMISSEquationsSetMaterialsCreateFinish
  
  !>Start the creation of materials for an equations set. \see OPENCMISS::CMISSEquationsSetMaterialsCreateFinish
  INTERFACE CMISSEquationsSetMaterialsCreateStart
    MODULE PROCEDURE CMISSEquationsSetMaterialsCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetMaterialsCreateStartObj
  END INTERFACE !CMISSEquationsSetMaterialsCreateStart
  
  !>Destroy the materials for an equations set.
  INTERFACE CMISSEquationsSetMaterialsDestroy
    MODULE PROCEDURE CMISSEquationsSetMaterialsDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetMaterialsDestroyObj
  END INTERFACE !CMISSEquationsSetMaterialsDestroy
  
  !>Returns the solution method for an equations set.
  INTERFACE CMISSEquationsSetSolutionMethodGet
    MODULE PROCEDURE CMISSEquationsSetSolutionMethodGetNumber
    MODULE PROCEDURE CMISSEquationsSetSolutionMethodGetObj
  END INTERFACE !CMISSEquationsSetSolutionMethodGet

  !>Sets/changes the solution method for an equations set.
  INTERFACE CMISSEquationsSetSolutionMethodSet
    MODULE PROCEDURE CMISSEquationsSetSolutionMethodSetNumber
    MODULE PROCEDURE CMISSEquationsSetSolutionMethodSetObj
  END INTERFACE !CMISSEquationsSetSolutionMethodSet
  
  !>Finish the creation of a source for an equations set. \see OPENCMISS::CMISSEquationsSetSourceCreateStart
  INTERFACE CMISSEquationsSetSourceCreateFinish
    MODULE PROCEDURE CMISSEquationsSetSourceCreateFinishNumber
    MODULE PROCEDURE CMISSEquationsSetSourceCreateFinishObj
  END INTERFACE !CMISSEquationsSetSourceCreateFinish
  
  !>Start the creation of a source for an equations set. \see OPENCMISS::CMISSEquationsSetSourceCreateFinish
  INTERFACE CMISSEquationsSetSourceCreateStart
    MODULE PROCEDURE CMISSEquationsSetSourceCreateStartNumber
    MODULE PROCEDURE CMISSEquationsSetSourceCreateStartObj
  END INTERFACE !CMISSEquationsSetSourceCreateStart
  
  !>Destroy the source for an equations set.
  INTERFACE CMISSEquationsSetSourceDestroy
    MODULE PROCEDURE CMISSEquationsSetSourceDestroyNumber
    MODULE PROCEDURE CMISSEquationsSetSourceDestroyObj
  END INTERFACE !CMISSEquationsSetSourceDestroy
  
  !>Returns the equations set specification i.e., equations set class, type and subtype for an equations set.
  INTERFACE CMISSEquationsSetSpecificationGet
    MODULE PROCEDURE CMISSEquationsSetSpecificationGetNumber
    MODULE PROCEDURE CMISSEquationsSetSpecificationGetObj
  END INTERFACE !CMISSEquationsSetSpecificationGet
  
   !>Sets/changes the equations set specification i.e., equations set class, type and subtype for an equations set.
  INTERFACE CMISSEquationsSetSpecificationSet
    MODULE PROCEDURE CMISSEquationsSetSpecificationSetNumber
    MODULE PROCEDURE CMISSEquationsSetSpecificationSetObj
  END INTERFACE !CMISSEquationsSetSpecificationSet
  
  PUBLIC CMISSEquationsSetAnalyticCreateFinish,CMISSEquationsSetAnalyticCreateStart
  
  PUBLIC CMISSEquationsSetAnalyticDestroy
  
  PUBLIC CMISSEquationsSetBoundaryConditionsAnalytic
  
  PUBLIC CMISSEquationsSetBoundaryConditionsCreateFinish,CMISSEquationsSetBoundaryConditionsCreateStart
  
  PUBLIC CMISSEquationsSetBoundaryConditionsDestroy
  
  PUBLIC CMISSEquationsSetCreateFinish,CMISSEquationsSetCreateStart
  
  PUBLIC CMISSEquationsSetDestroy
  
  PUBLIC CMISSEquationsSetDependentCreateFinish,CMISSEquationsSetDependentCreateStart
  
  PUBLIC CMISSEquationsSetDependentDestroy
  
  PUBLIC CMISSEquationsSetEquationsCreateFinish,CMISSEquationsSetEquationsCreateStart
  
  PUBLIC CMISSEquationsSetEquationsDestroy
  
  PUBLIC CMISSEquationsSetMaterialsCreateFinish,CMISSEquationsSetMaterialsCreateStart
  
  PUBLIC CMISSEquationsSetMaterialsDestroy

  PUBLIC CMISSEquationsSetSolutionMethodGet,CMISSEquationsSetSolutionMethodSet
  
  PUBLIC CMISSEquationsSetSourceCreateFinish,CMISSEquationsSetSourceCreateStart
  
  PUBLIC CMISSEquationsSetSourceDestroy

  PUBLIC CMISSEquationsSetSpecificationGet,CMISSEquationsSetSpecificationSet
  
  
!!==================================================================================================================================
!!
!! FIELD_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !> \addtogroup OPENCMISS_FieldConstants OPENCMISS::Field::Constants
  !> \brief Field constants.
  !>@{  
  !> \addtogroup OPENCMISS_FieldDependentTypes OPENCMISS::Field::DependentTypes
  !> \brief Depedent field parameter types.
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldIndependentType = FIELD_INDEPENDENT_TYPE !<Independent field type \see OPENCMISS_FieldDependentTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDependentType = FIELD_DEPENDENT_TYPE !<Dependent field type \see OPENCMISS_FieldDependentTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldDimensionTypes OPENCMISS::Field::DimensionTypes
  !> \brief Field dimension parameter types.
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldScalarDimensionType = FIELD_SCALAR_DIMENSION_TYPE !<Scalar field \see OPENCMISS_FieldDimensionTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldVectorDimensionType = FIELD_VECTOR_DIMENSION_TYPE !<Vector field \see OPENCMISS_FieldDimensionTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldTensorDimensionType = FIELD_TENSOR_DIMENSION_TYPE !<Tensor field \see OPENCMISS_FieldDimensionTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldTypes OPENCMISS::Field::Types
  !> \brief Field type parameters.
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldGeometricType = FIELD_GEOMETRIC_TYPE !<Geometric field \see OPENCMISS_FieldTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldFibreType = FIELD_FIBRE_TYPE !<Fibre field \see OPENCMISS_FieldTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldGeneralType = FIELD_GENERAL_TYPE !<General field \see OPENCMISS_FieldTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldMaterialType = FIELD_MATERIAL_TYPE !<Material field \see OPENCMISS_FieldTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldInterpolationTypes OPENCMISS::Field::InterpolationTypes
  !> \brief Field interpolation parameters.
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldConstantInterpolation = FIELD_CONSTANT_INTERPOLATION !<Constant interpolation. One parameter for the field \see OPENCMISS_FieldInterpolationTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldElementBasedInterpolation = FIELD_ELEMENT_BASED_INTERPOLATION !<Element based interpolation. Parameters are different in each element \see OPENCMISS_FieldInterpolationTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldNodeBasedInterpolation = FIELD_NODE_BASED_INTERPOLATION !<Node based interpolation. Parameters are nodal based and a basis function is used \see OPENCMISS_FieldInterpolationTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldGridPointBasedInterpolation = FIELD_GRID_POINT_BASED_INTERPOLATION !<Grid point based interpolation. Parameters are different at each grid point \see OPENCMISS_FieldInterpolationTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldGaussPointBasedInterpolation = FIELD_GAUSS_POINT_BASED_INTERPOLATION !<Gauss point based interpolation. Parameters are different at each Gauss point \see OPENCMISS_FieldInterpolationTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldVariableTypes OPENCMISS::Field::VariableTypes
  !> \brief Field variable type parameters.
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldUVariableType = FIELD_U_VARIABLE_TYPE !<Standard variable type i.e., u \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDelUDelNVariableType = FIELD_DELUDELN_VARIABLE_TYPE !<Normal derivative variable type i.e., du/dn \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDelUDelTVariableType = FIELD_DELUDELT_VARIABLE_TYPE !<First time derivative variable type i.e., du/dt \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDel2UDelT2VariableType = FIELD_DEL2UDELT2_VARIABLE_TYPE !<Second type derivative variable type i.e., d^2u/dt^2 \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldVVariableType = FIELD_V_VARIABLE_TYPE !<Second standard variable type i.e., v \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDelVDelNVariableType = FIELD_DELVDELN_VARIABLE_TYPE !<Second normal variable type i.e., dv/dn \see OPENCMISS_FieldVariableTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldDataTypes OPENCMISS::Field::DataTypes
  !> \brief Field data types
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldIntgType = FIELD_INTG_TYPE !<Integer field data type \see OPENCMISS_FieldDataTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldSPType = FIELD_SP_TYPE !<Single precision real field data type \see OPENCMISS_FieldDataTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDPType = FIELD_DP_TYPE !<Double precision real field data type \see OPENCMISS_FieldDataTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldLType = FIELD_L_TYPE !<Logical field data type \see OPENCMISS_FieldDataTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldDOFOrderTypes OPENCMISS::Field::DOFOrderTypes
  !> \brief Field DOF order types
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldSeparatedComponentDOFOrder = FIELD_SEPARATED_COMPONENT_DOF_ORDER !<Field variable component dofs are not contiguous \see OPENCMISS_FieldDOFOrderTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldContiguousComponentDOFOrder = FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER !<Field variable component dofs are contiguous \see OPENCMISS_FieldDOFOrderTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldParameterSetTypes OPENCMISS::Field::ParameterSetTypes
  !> \brief Field parameter set type parameters
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldValuesSetType = FIELD_VALUES_SET_TYPE !<The parameter set corresponding to the field values (at time T+DT for dynamic problems) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSInitialValuesSetType = FIELD_INITIAL_VALUES_SET_TYPE !<The parameter set corresponding to the field initial values \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldIncrementalValuesSetType = FIELD_INCREMENTAL_VALUES_SET_TYPE !<The parameter set corresponding to the field incremental values \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldAnalyticValuesSetType = FIELD_ANALYTIC_VALUES_SET_TYPE !<The parameter set corresponding to the analytic field values \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSPreviousValuesSetType = FIELD_PREVIOUS_VALUES_SET_TYPE !<The parameter set corresponding to the previous field values (at time T) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSMeanPredictedDisplacementSetType = FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE !<The parameter set corresponding to the mean predicited avalues (at time T+DT) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldVelocityValuesSetType = FIELD_VELOCITY_VALUES_SET_TYPE !<The parameter set corresponding to the velocity values (at time T+DT) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldInitialVelocitySetType = FIELD_INITIAL_VELOCITY_SET_TYPE !<The parameter set corresponding to the initial velocity values for dynamic problems. This is also the previous velocity values \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldPreviousVelocitySetType = FIELD_PREVIOUS_VELOCITY_SET_TYPE !<The parameter set corresponding to the previous velocity values (at time T). This is also the initial velocity values for dynamic problems. \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldMeanPredictedVelocitySetType = FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE !<The parameter set corresponding to the mean predicited velocity values (at time T+DT) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldAccelerationValuesSetType = FIELD_ACCELERATION_VALUES_SET_TYPE !<The parameter set corresponding to the acceleration values (at time T+DT) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSInitialAccelerationSetType = FIELD_INITIAL_ACCELERATION_SET_TYPE !<The parameter set corresponding to the initial acceleration values for dynamic problems. This is also the previous accelearation values \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldPreviousAccelerationSetType = FIELD_PREVIOUS_ACCELERATION_SET_TYPE !<The parameter set corresponding to the previous acceleration values (at time T).This is also the initial acceleration values for dynamic problems. \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSMeanPredictedAccelerationSetType = FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE !<The parameter set corresponding to the mean predicited acceleration values (at time T+DT) \see OPENCMISS_FieldParameterSetTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_FieldScalingTypes OPENCMISS::Field::ScalingTypes
  !> \brief Field scaling type parameters
  !> \see OPENCMISS::Field,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSFieldNoScaling = FIELD_NO_SCALING !<The field is not scaled \see OPENCMISS_FieldScalingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldUnitScaling = FIELD_UNIT_SCALING !<The field has unit scaling \see OPENCMISS_FieldScalingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldArcLengthScaling = FIELD_ARC_LENGTH_SCALING !<The field has arc length scaling \see OPENCMISS_FieldScalingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldArithmeticMeanScaling = FIELD_ARITHMETIC_MEAN_SCALING !<The field has arithmetic mean of the arc length scaling \see OPENCMISS_FieldScalingTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldHarmonicMeanScaling = FIELD_HARMONIC_MEAN_SCALING !<The field has geometric mean of the arc length scaling \see OPENCMISS_FieldScalingTypes,OPENCMISS
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Returns the interpolation type for a field variable component.
  INTERFACE CMISSFieldComponentInterpolationGet
    MODULE PROCEDURE CMISSFieldComponentInterpolationGetNumber
    MODULE PROCEDURE CMISSFieldComponentInterpolationGetObj
  END INTERFACE !CMISSFieldComponentInterpolationGet

  !>Sets/changes the interpolation type for a field variable component.
  INTERFACE CMISSFieldComponentInterpolationSet
    MODULE PROCEDURE CMISSFieldComponentInterpolationSetNumber
    MODULE PROCEDURE CMISSFieldComponentInterpolationSetObj
  END INTERFACE !CMISSFieldComponentInterpolationSet

  !>Returns the label for a field variable component.
  INTERFACE CMISSFieldComponentLabelGet
    MODULE PROCEDURE CMISSFieldComponentLabelGetCNumber
    MODULE PROCEDURE CMISSFieldComponentLabelGetCObj
    MODULE PROCEDURE CMISSFieldComponentLabelGetVSNumber
    MODULE PROCEDURE CMISSFieldComponentLabelGetVSObj
  END INTERFACE !CMISSFieldComponentLabelGet
  
  !>Sets/changes the label for a field variable component.
  INTERFACE CMISSFieldComponentLabelSet
    MODULE PROCEDURE CMISSFieldComponentLabelSetCNumber
    MODULE PROCEDURE CMISSFieldComponentLabelSetCObj
    MODULE PROCEDURE CMISSFieldComponentLabelSetVSNumber
    MODULE PROCEDURE CMISSFieldComponentLabelSetVSObj
  END INTERFACE !CMISSFieldComponentLabelSet

  !>Returns the mesh component number for a field variable component.
  INTERFACE CMISSFieldComponentMeshComponentGet
    MODULE PROCEDURE CMISSFieldComponentMeshComponentGetNumber
    MODULE PROCEDURE CMISSFieldComponentMeshComponentGetObj
  END INTERFACE !CMISSFieldComponentMeshComponentGet

  !>Sets/changes the mesh component number for a field variable component.
  INTERFACE CMISSFieldComponentMeshComponentSet
    MODULE PROCEDURE CMISSFieldComponentMeshComponentSetNumber
    MODULE PROCEDURE CMISSFieldComponentMeshComponentSetObj
  END INTERFACE !CMISSFieldComponentMeshComponentSet  

  !>Initialises the values of a parameter set of a field variable component to a constant value.
  INTERFACE CMISSFieldComponentValuesInitialise
    MODULE PROCEDURE CMISSFieldComponentValuesInitialiseNumber
    MODULE PROCEDURE CMISSFieldComponentValuesInitialiseObj
  END INTERFACE !CMISSFieldComponentValuesInitialise

  !>Returns the data type for a field variable.
  INTERFACE CMISSFieldDataTypeGet
    MODULE PROCEDURE CMISSFieldDataTypeGetNumber
    MODULE PROCEDURE CMISSFieldDataTypeGetObj
  END INTERFACE !CMISSFieldDataTypeGet

  !>Sets/changes the data type for a field variable.
  INTERFACE CMISSFieldDataTypeSet
    MODULE PROCEDURE CMISSFieldDataTypeSetNumber
    MODULE PROCEDURE CMISSFieldDataTypeSetObj
  END INTERFACE !CMISSFieldDataTypeSet

  !>Returns the DOF order type for a field variable.
  INTERFACE CMISSFieldDOFOrderTypeGet
    MODULE PROCEDURE CMISSFieldDOFOrderTypeGetNumber
    MODULE PROCEDURE CMISSFieldDOFOrderTypeGetObj
  END INTERFACE !CMISSFieldDOFOrderTypeGet

  !>Sets/changes the DOF order type for a field variable. Note: for contiguous coponent DOF ordering all the components of the field variable must have the same interpolation type.
  INTERFACE CMISSFieldDOFOrderTypeSet
    MODULE PROCEDURE CMISSFieldDOFOrderTypeSetNumber
    MODULE PROCEDURE CMISSFieldDOFOrderTypeSetObj
  END INTERFACE !CMISSFieldDOFOrderTypeSet

  !>Finishes the creation of a field. \see OPENCMISS::CMISSFieldCreateStart
  INTERFACE CMISSFieldCreateFinish 
    MODULE PROCEDURE CMISSFieldCreateFinishNumber
    MODULE PROCEDURE CMISSFieldCreateFinishObj
  END INTERFACE !CMISSFieldCreateFinish

  !>Starts the creation of a field. \see OPENCMISS::CMISSFieldCreateFinish
  INTERFACE CMISSFieldCreateStart
    MODULE PROCEDURE CMISSFieldCreateStartNumber
    MODULE PROCEDURE CMISSFieldCreateStartObj
  END INTERFACE !CMISSFieldCreateStart

  !>Returns the dependent type for a field.
  INTERFACE CMISSFieldDependentTypeGet
    MODULE PROCEDURE CMISSFieldDependentTypeGetNumber
    MODULE PROCEDURE CMISSFieldDependentTypeGetObj
  END INTERFACE !CMISSFieldDependentTypeGet

  !>Sets/changes the dependent type for a field.
  INTERFACE CMISSFieldDependentTypeSet
    MODULE PROCEDURE CMISSFieldDependentTypeSetNumber
    MODULE PROCEDURE CMISSFieldDependentTypeSetObj
  END INTERFACE !CMISSFieldDependentTypeSet

  !>Destroys a field.
  INTERFACE CMISSFieldDestroy
    MODULE PROCEDURE CMISSFieldDestroyNumber
    MODULE PROCEDURE CMISSFieldDestroyObj
  END INTERFACE !CMISSFieldDestroy

  !>Returns the field dimension for a field variable.
  INTERFACE CMISSFieldDimensionGet
    MODULE PROCEDURE CMISSFieldDimensionGetNumber
    MODULE PROCEDURE CMISSFieldDimensionGetObj
  END INTERFACE !CMISSFieldDimensionGet

  !>Sets/changes the field dimension for a field variable.
  INTERFACE CMISSFieldDimensionSet
    MODULE PROCEDURE CMISSFieldDimensionSetNumber
    MODULE PROCEDURE CMISSFieldDimensionSetObj
  END INTERFACE !CMISSFieldDimensionSet

  !>Returns the geometric field for a field.
  INTERFACE CMISSFieldGeometricFieldGet
    MODULE PROCEDURE CMISSFieldGeometricFieldGetNumber
    MODULE PROCEDURE CMISSFieldGeometricFieldGetObj
  END INTERFACE !CMISSFieldGeometricFieldGet

  !>Sets/changes the geometric field for a field. 
  INTERFACE CMISSFieldGeometricFieldSet
    MODULE PROCEDURE CMISSFieldGeometricFieldSetNumber
    MODULE PROCEDURE CMISSFieldGeometricFieldSetObj
  END INTERFACE !CMISSFieldGeometricFieldSet

 !>Returns the label for a field.
  INTERFACE CMISSFieldLabelGet
    MODULE PROCEDURE CMISSFieldLabelGetCNumber
    MODULE PROCEDURE CMISSFieldLabelGetCObj
    MODULE PROCEDURE CMISSFieldLabelGetVSNumber
    MODULE PROCEDURE CMISSFieldLabelGetVSObj
  END INTERFACE !CMISSFieldLabelGet
  
  !>Sets/changes the label for a field.
  INTERFACE CMISSFieldLabelSet
    MODULE PROCEDURE CMISSFieldLabelSetCNumber
    MODULE PROCEDURE CMISSFieldLabelSetCObj
    MODULE PROCEDURE CMISSFieldLabelSetVSNumber
    MODULE PROCEDURE CMISSFieldLabelSetVSObj
  END INTERFACE !CMISSFieldLabelSet

  !>Returns the mesh decomposition for a field. 
  INTERFACE CMISSFieldMeshDecompositionGet
    MODULE PROCEDURE CMISSFieldMeshDecompositionGetNumber
    MODULE PROCEDURE CMISSFieldMeshDecompositionGetObj
  END INTERFACE !CMISSFieldMeshDecompositionGet

  !>Sets/changes the mesh decomposition for a field. \todo remove when fields take decomposition argument on creation???
  INTERFACE CMISSFieldMeshDecompositionSet
    MODULE PROCEDURE CMISSFieldMeshDecompositionSetNumber
    MODULE PROCEDURE CMISSFieldMeshDecompositionSetObj
  END INTERFACE !CMISSFieldMeshDecompositionSet

  !>Returns the number of field components for a field variable.
  INTERFACE CMISSFieldNumberOfComponentsGet 
    MODULE PROCEDURE CMISSFieldNumberOfComponentsGetNumber
    MODULE PROCEDURE CMISSFieldNumberOfComponentsGetObj
  END INTERFACE !CMISSFieldNumberOfComponentsGet

  !>Sets/changes the number of field components for a field variable.
  INTERFACE CMISSFieldNumberOfComponentsSet 
    MODULE PROCEDURE CMISSFieldNumberOfComponentsSetNumber
    MODULE PROCEDURE CMISSFieldNumberOfComponentsSetObj
  END INTERFACE !CMISSFieldNumberOfComponentsSet

  !>Returns the number of field variables for a field.
  INTERFACE CMISSFieldNumberOfVariablesGet 
    MODULE PROCEDURE CMISSFieldNumberOfVariablesGetNumber
    MODULE PROCEDURE CMISSFieldNumberOfVariablesGetObj
  END INTERFACE !CMISSFieldNumberOfVariablesGet

  !>Sets/changes the number of field variables for a field.
  INTERFACE CMISSFieldNumberOfVariablesSet 
    MODULE PROCEDURE CMISSFieldNumberOfVariablesSetNumber
    MODULE PROCEDURE CMISSFieldNumberOfVariablesSetObj
  END INTERFACE !CMISSFieldNumberOfVariablesSet

  !>Adds the given value to the given parameter set for the constant of the field variable component.
  INTERFACE CMISSFieldParameterSetAddConstant
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantSPObj
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantDPObj
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantLNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddConstantLObj
  END INTERFACE !CMISSFieldParameterSetAddConstant

  !>Adds the given value to the given parameter set for a particular user element of the field variable component.
  INTERFACE CMISSFieldParameterSetAddElement
    MODULE PROCEDURE CMISSFieldParameterSetAddElementIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddElementIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetAddElementSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddElementSPObj
    MODULE PROCEDURE CMISSFieldParameterSetAddElementDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddElementDPObj
    MODULE PROCEDURE CMISSFieldParameterSetAddElementLNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddElementLObj
  END INTERFACE !CMISSFieldParameterSetAddElement

  !>Adds the given value to the given parameter set for a particular user node of the field variable component.
  INTERFACE CMISSFieldParameterSetAddNode
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeSPObj
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeDPObj
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeLNumber
    MODULE PROCEDURE CMISSFieldParameterSetAddNodeLObj
  END INTERFACE !CMISSFieldParameterSetAddNode

  !>Creates a new parameter set of type set type for a field variable.
  INTERFACE CMISSFieldParameterSetCreate
    MODULE PROCEDURE CMISSFieldParameterSetCreateNumber
    MODULE PROCEDURE CMISSFieldParameterSetCreateObj
  END INTERFACE !CMISSFieldParameterSetCreate
  
  !>Destroy a parameter set of type set type for a field variable.
  INTERFACE CMISSFieldParameterSetDestroy
    MODULE PROCEDURE CMISSFieldParameterSetDestroyNumber
    MODULE PROCEDURE CMISSFieldParameterSetDestroyObj
  END INTERFACE !CMISSFieldParameterSetCreate
  
  !>Returns a pointer to the specified field parameter set local data array. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  INTERFACE CMISSFieldParameterSetDataGet
    MODULE PROCEDURE CMISSFieldParameterSetDataGetIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataGetIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetDataGetSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataGetSPObj
    MODULE PROCEDURE CMISSFieldParameterSetDataGetDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataGetDPObj
    MODULE PROCEDURE CMISSFieldParameterSetDataGetLNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataGetLObj
  END INTERFACE !CMISSFieldParameterSetDataGet
  
  !>Restores the specified field variable parameter set local array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call.
  INTERFACE CMISSFieldParameterSetDataRestore
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreSPObj
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreDPObj
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreLNumber
    MODULE PROCEDURE CMISSFieldParameterSetDataRestoreLObj
  END INTERFACE !CMISSFieldParameterSetDataRestore
  
  !>Updates the given parameter set with the given value for the constant of a field variable component.
  INTERFACE CMISSFieldParameterSetUpdateConstant
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantSPObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantDPObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantLNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateConstantLObj
  END INTERFACE !CMISSFieldParameterSetUpdateConstant

  !>Updates the given parameter set with the given value for a particular user element of a field variable component.
  INTERFACE CMISSFieldParameterSetUpdateElement
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementSPObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementDPObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementLNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateElementLObj
  END INTERFACE !CMISSFieldParameterSetUpdateElement

  !>Finishes the parameter set update for a field variable. \see OPENCMISS::CMISSFieldParameterSetUpdateStart
  INTERFACE CMISSFieldParameterSetUpdateFinish
    MODULE PROCEDURE CMISSFieldParameterSetUpdateFinishNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateFinishObj
  END INTERFACE !CMISSFieldParameterSetUpdateFinish

  !>Updates the given parameter set with the given value for a particular user node of a field variable component.
  INTERFACE CMISSFieldParameterSetUpdateNode
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeIntgNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeIntgObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeSPNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeSPObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeDPNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeDPObj
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeLNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateNodeLObj
  END INTERFACE !CMISSFieldParameterSetUpdateNode

  !>Starts the parameter set update for a field variable. \see OPENCMISS::CMISSFieldParameterSetUpdateFinish
  INTERFACE CMISSFieldParameterSetUpdateStart
    MODULE PROCEDURE CMISSFieldParameterSetUpdateStartNumber
    MODULE PROCEDURE CMISSFieldParameterSetUpdateStartObj
  END INTERFACE !CMISSFieldParameterSetUpdateStart

  !>Returns the scaling type for a field.  
  INTERFACE CMISSFieldScalingTypeGet
    MODULE PROCEDURE CMISSFieldScalingTypeGetNumber
    MODULE PROCEDURE CMISSFieldScalingTypeGetObj
  END INTERFACE !CMISSFieldScalingTypeGet
    
  !>Sets/changes the scaling type for a field.  
  INTERFACE CMISSFieldScalingTypeSet
    MODULE PROCEDURE CMISSFieldScalingTypeSetNumber
    MODULE PROCEDURE CMISSFieldScalingTypeSetObj
  END INTERFACE !CMISSFieldScalingTypeSet
    
  !>Returns the type for a field. 
  INTERFACE CMISSFieldTypeGet
    MODULE PROCEDURE CMISSFieldTypeGetNumber
    MODULE PROCEDURE CMISSFieldTypeGetObj
  END INTERFACE !CMISSFieldTypeGet
    
  !>Sets/changes the type for a field. 
  INTERFACE CMISSFieldTypeSet
    MODULE PROCEDURE CMISSFieldTypeSetNumber
    MODULE PROCEDURE CMISSFieldTypeSetObj
  END INTERFACE !CMISSFieldTypeSet
    
  !>Returns the label for a field variable. 
  INTERFACE CMISSFieldVariableLabelGet
    MODULE PROCEDURE CMISSFieldVariableLabelGetCNumber
    MODULE PROCEDURE CMISSFieldVariableLabelGetCObj
    MODULE PROCEDURE CMISSFieldVariableLabelGetVSNumber
    MODULE PROCEDURE CMISSFieldVariableLabelGetVSObj
  END INTERFACE !CMISSFieldVariableLabelGet
    
  !>Sets/changes the label for a field variable. 
  INTERFACE CMISSFieldVariableLabelSet
    MODULE PROCEDURE CMISSFieldVariableLabelSetCNumber
    MODULE PROCEDURE CMISSFieldVariableLabelSetCObj
    MODULE PROCEDURE CMISSFieldVariableLabelSetVSNumber
    MODULE PROCEDURE CMISSFieldVariableLabelSetVSObj
  END INTERFACE !CMISSFieldVariableLabelSet
    
  !>Returns the field variable types for a field. 
  INTERFACE CMISSFieldVariableTypesGet
    MODULE PROCEDURE CMISSFieldVariableTypesGetNumber
    MODULE PROCEDURE CMISSFieldVariableTypesGetObj
  END INTERFACE !CMISSFieldVariableTypesGet
    
  !>Sets/changes the field variable types for a field. 
  INTERFACE CMISSFieldVariableTypesSet
    MODULE PROCEDURE CMISSFieldVariableTypesSetNumber
    MODULE PROCEDURE CMISSFieldVariableTypesSetObj
  END INTERFACE !CMISSFieldVariableTypesSet
    
  PUBLIC CMISSFieldDependentType,CMISSFieldIndependentType

  PUBLIC CMISSFieldScalarDimensionType,CMISSFieldVectorDimensionType,CMISSFieldTensorDimensionType

  PUBLIC CMISSFieldGeometricType,CMISSFieldFibreType,CMISSFieldGeneralType,CMISSFieldMaterialType

  PUBLIC CMISSFieldConstantInterpolation,CMISSFieldElementBasedInterpolation,CMISSFieldNodeBasedInterpolation, &
    & CMISSFieldGridPointBasedInterpolation,CMISSFieldGaussPointBasedInterpolation

  PUBLIC CMISSFieldUVariableType,CMISSFieldDelUDelNVariableType,CMISSFieldDelUDelTVariableType,CMISSFieldDel2UDelT2VariableType, &
    & CMISSFieldVVariableType,CMISSFieldDelVDelNVariableType

  PUBLIC CMISSFieldIntgType,CMISSFieldSPType,CMISSFieldDPType,CMISSFieldLType

  PUBLIC CMISSFieldSeparatedComponentDOFOrder,CMISSFieldContiguousComponentDOFOrder

  PUBLIC CMISSFieldValuesSetType,CMISSInitialValuesSetType,CMISSFieldIncrementalValuesSetType,CMISSFieldAnalyticValuesSetType, &
    & CMISSPreviousValuesSetType,CMISSMeanPredictedDisplacementSetType,CMISSFieldVelocityValuesSetType, &
    & CMISSFieldInitialVelocitySetType,CMISSFieldPreviousVelocitySetType,CMISSFieldMeanPredictedVelocitySetType, &
    & CMISSFieldAccelerationValuesSetType,CMISSInitialAccelerationSetType,CMISSFieldPreviousAccelerationSetType, &
    & CMISSMeanPredictedAccelerationSetType

  PUBLIC CMISSFieldNoScaling,CMISSFieldUnitScaling,CMISSFieldArcLengthScaling,CMISSFieldArithmeticMeanScaling, &
    & CMISSFieldHarmonicMeanScaling

  PUBLIC CMISSFieldComponentInterpolationGet,CMISSFieldComponentInterpolationSet

  PUBLIC CMISSFieldComponentLabelGet,CMISSFieldComponentLabelSet

  PUBLIC CMISSFieldComponentMeshComponentGet,CMISSFieldComponentMeshComponentSet

  PUBLIC CMISSFieldComponentValuesInitialise

  PUBLIC CMISSFieldDataTypeGet,CMISSFieldDataTypeSet

  PUBLIC CMISSFieldDOFOrderTypeGet,CMISSFieldDOFOrderTypeSet

  PUBLIC CMISSFieldCreateFinish,CMISSFieldCreateStart

  PUBLIC CMISSFieldDependentTypeGet,CMISSFieldDependentTypeSet

  PUBLIC CMISSFieldDestroy

  PUBLIC CMISSFieldDimensionGet,CMISSFieldDimensionSet

  PUBLIC CMISSFieldGeometricFieldGet,CMISSFieldGeometricFieldSet

  PUBLIC CMISSFieldLabelGet,CMISSFieldLabelSet

  PUBLIC CMISSFieldMeshDecompositionGet,CMISSFieldMeshDecompositionSet

  PUBLIC CMISSFieldNumberOfComponentsGet,CMISSFieldNumberOfComponentsSet

  PUBLIC CMISSFieldNumberOfVariablesGet,CMISSFieldNumberOfVariablesSet

  PUBLIC CMISSFieldParameterSetAddConstant,CMISSFieldParameterSetAddElement,CMISSFieldParameterSetAddNode

  PUBLIC CMISSFieldParameterSetCreate

  PUBLIC CMISSFieldParameterSetDestroy

  PUBLIC CMISSFieldParameterSetDataGet,CMISSFieldParameterSetDataRestore

  PUBLIC CMISSFieldParameterSetUpdateConstant,CMISSFieldParameterSetUpdateElement,CMISSFieldParameterSetUpdateNode

  PUBLIC CMISSFieldParameterSetUpdateFinish,CMISSFieldParameterSetUpdateStart

  PUBLIC CMISSFieldScalingTypeGet,CMISSFieldScalingTypeSet

  PUBLIC CMISSFieldTypeGet,CMISSFieldTypeSet

  PUBLIC CMISSFieldVariableLabelGet,CMISSFieldVariableLabelSet

  PUBLIC CMISSFieldVariableTypesGet,CMISSFieldVariableTypesSet
 
!!==================================================================================================================================
!!
!! FIELD_IO_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE CMISSFieldIOElementsExport
    MODULE PROCEDURE CMISSFieldIOElementsExportObj
  END INTERFACE !CMISSFieldIOElementsExport

  INTERFACE CMISSFieldIONodesExport
    MODULE PROCEDURE CMISSFieldIONodesExportObj
  END INTERFACE !CMISSFieldIONodesExport

  PUBLIC CMISSFieldIOElementsExport,CMISSFieldIONodesExport
  
!!==================================================================================================================================
!!
!! PROBLEM_CONSTANTS_ROUTINES
!!
!!==================================================================================================================================

  !Module parameters
 
  !> \addtogroup OPENCMISS_ProblemConstants OPENCMISS::Problem::Constants
  !> \brief Problem constants.
  !>@{  
  !> \addtogroup OPENCMISS_ProblemClasses OPENCMISS::Problem::Classes
  !> \brief Problem classes.
  !> \see OPENCMISS::Problem,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSProblemNoClass = PROBLEM_NO_CLASS !<No problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemElasticityClass = PROBLEM_ELASTICITY_CLASS !<Elasticity problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemFluidMechanicsClass = PROBLEM_FLUID_MECHANICS_CLASS !<Fluid mechanics problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemElectromagneticsClass = PROBLEM_ELECTROMAGNETICS_CLASS !<Electromagnetics problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemClassicalFieldClass = PROBLEM_CLASSICAL_FIELD_CLASS !<Classical field problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemBioelectricsClass = PROBLEM_BIOELECTRICS_CLASS !<Bioelectrics problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemModalClass = PROBLEM_MODAL_CLASS !<Modal problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemFittingClass = PROBLEM_FITTING_CLASS !<Fitting problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemOptimisationClass = PROBLEM_OPTIMISATION_CLASS !<Optimisation problem class \see OPENCMISS_ProblemClasses,OPENCMISS 
  !>@}
  !> \addtogroup OPENCMISS_ProblemTypes OPENCMISS::Problem::Types
  !> \brief Problem Types.
  !> \see OPENCMISS::Problem,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSProblemNoType = PROBLEM_NO_TYPE !<No problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemLinearElasticityType = PROBLEM_LINEAR_ELASTICITY_TYPE !<Linear elasticity problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemFiniteElasticityType = PROBLEM_FINITE_ELASTICITY_TYPE !<Finite elasticity problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemStokesEquationType = PROBLEM_STOKES_EQUATION_TYPE !<Stokes equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemNavierStokesType = PROBLEM_NAVIER_STOKES_TYPE !<Navier-Stokes problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemDarcyEquationType = PROBLEM_DARCY_EQUATION_TYPE !<Darcy equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemElectrostaticType = PROBLEM_ELECTROSTATIC_TYPE !<Electrostatic problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemMagnetostaticType = PROBLEM_MAGNETOSTATIC_TYPE !<Magnetostatic problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemMaxwellsEquationsType = PROBLEM_MAXWELLS_EQUATIONS_TYPE !<Maxwell's equations problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemLaplaceEquationType = PROBLEM_LAPLACE_EQUATION_TYPE !<Laplace problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemPoissonEquationType = PROBLEM_POISSON_EQUATION_TYPE !<Poisson problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemHelmholtzEquationType = PROBLEM_HELMHOLTZ_EQUATION_TYPE !<Helmholtz problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemWaveEquationType = PROBLEM_WAVE_EQUATION_TYPE !<Wave equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemDiffusionEquationType = PROBLEM_DIFFUSION_EQUATION_TYPE !<Diffusion equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemAdvectionDiffusionEquationType = PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE !<Advection-Diffusion equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemReactionDiffusionEquationType = PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE !<Reaction-Diffusion equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemBiharmonicEquationType = PROBLEM_BIHARMONIC_EQUATION_TYPE !<Bi-harmonic equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemMonodomainEquationType = PROBLEM_MONODOMAIN_EQUATION_TYPE !<Monodomain equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemBidomainEquationType = PROBLEM_BIDOMAIN_EQUATION_TYPE !<Bidomain equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
  INTEGER(INTG), PARAMETER :: CMISSProblemLinearElasticModalType = PROBLEM_LINEAR_ELASTIC_MODAL_TYPE !<Linear elastic modal problem type \see OPENCMISS_ProblemTypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_ProblemSubTypes OPENCMISS::Problem::Subtypes
  !> \brief Problem Subtypes.
  !> \see OPENCMISS::Problem,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSProblemNoSubtype = PROBLEM_NO_SUBTYPE !<No problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemStaticStokesSubtype = PROBLEM_STATIC_STOKES_SUBTYPE !<Static Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemLaplaceStokesSubtype = PROBLEM_LAPLACE_STOKES_SUBTYPE !<Laplace type Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemTransientStokesSubtype = PROBLEM_TRANSIENT_STOKES_SUBTYPE !<Transient Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemOptimisedSotkesSubtype = PROBLEM_OPTIMISED_STOKES_SUBTYPE !<Optimised Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemStaticNavierStokesSubtype = PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE !<Static Navier-Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemLaplaceNavierStokesSubtype = PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE !<Laplace type Navier-Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemTransientNavierStokesSubtype = PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE !<Transient Navier-Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemOptimisedNavierStokesSubtype = PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE !<Optimised Navier-Stokes problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemStandardDarcySubtype = PROBLEM_STANDARD_DARCY_SUBTYPE !<Standard Darcy problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemStandardLaplaceSubtype = PROBLEM_STANDARD_LAPLACE_SUBTYPE !<Standard Laplace problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemGeneralisedLaplaceSubtype = PROBLEM_GENERALISED_LAPLACE_SUBTYPE !<Generalised Laplace problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemLinearSourcePoissonSubtype = PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE !<Linear source Poisson problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemNonlinearSourcePoissonSubtype = PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE !<Nonlinear source Poisson problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemNoSourceHelmholtzSubtype = PROBLEM_NO_SOURCE_HELMHOLTZ_SUBTYPE !<No source Helmholtz problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemNoSourceDiffusionSubtype = PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE !<No source Diffusion problem subtype \see OPENCMISS_ProblemSubtypes,OPENCMISS
  !>@}
  !> \addtogroup OPENCMISS_ProblemControlLoopTypes OPENCMISS::Problem::ControlLoopTypes
  !> \brief Problem control loop type parameters
  !> \see OPENCMISS::Problem,OPENCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISSProblemControlSimpleType = PROBLEM_CONTROL_SIMPLE_TYPE !<Simple, one iteration control loop. \see OPENCMISS_ProblemControlLoopTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemControlFixedLoopType = PROBLEM_CONTROL_FIXED_LOOP_TYPE !<Fixed iteration control loop. \see OPENCMISS_ProblemControlLoopTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemControlTimeLoopType = PROBLEM_CONTROL_TIME_LOOP_TYPE !<Time control loop. \see OPENCMISS_ProblemControlLoopTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSProblemControlWhileLoopType = PROBLEM_CONTROL_WHILE_LOOP_TYPE !<While control loop. \see OPENCMISS_ProblemControlLoopTypes,OPENCMISS
  !>@}
  !>@}

  PUBLIC CMISSProblemNoClass,CMISSProblemElasticityClass,CMISSProblemFluidMechanicsClass,CMISSProblemElectromagneticsClass, &
    & CMISSProblemClassicalFieldClass,CMISSProblemBioelectricsClass,CMISSProblemModalClass,CMISSProblemFittingClass, &
    & CMISSProblemOptimisationClass

  PUBLIC CMISSProblemNoType

  PUBLIC CMISSProblemLinearElasticityType,CMISSProblemFiniteElasticityType

  PUBLIC CMISSProblemStokesEquationType,CMISSProblemNavierStokesType,CMISSProblemDarcyEquationType

  PUBLIC CMISSProblemElectrostaticType,CMISSProblemMagnetostaticType,CMISSProblemMaxwellsEquationsType

  PUBLIC CMISSProblemLaplaceEquationType,CMISSProblemPoissonEquationType,CMISSProblemHelmholtzEquationType, &
    & CMISSProblemWaveEquationType,CMISSProblemDiffusionEquationType,CMISSProblemAdvectionDiffusionEquationType, &
    & CMISSProblemReactionDiffusionEquationType,CMISSProblemBiharmonicEquationType

  PUBLIC CMISSProblemMonodomainEquationType,CMISSProblemBidomainEquationType

  PUBLIC CMISSProblemLinearElasticModalType

  PUBLIC CMISSProblemNoSubtype

  PUBLIC CMISSProblemStaticStokesSubtype,CMISSProblemLaplaceStokesSubtype,CMISSProblemTransientStokesSubtype, &
    & CMISSProblemOptimisedSotkesSubtype

  PUBLIC CMISSProblemStaticNavierStokesSubtype,CMISSProblemLaplaceNavierStokesSubtype,CMISSProblemTransientNavierStokesSubtype, &
    & CMISSProblemOptimisedNavierStokesSubtype

  PUBLIC CMISSProblemStandardDarcySubtype

  PUBLIC CMISSProblemStandardLaplaceSubtype,CMISSProblemGeneralisedLaplaceSubtype

  PUBLIC CMISSProblemLinearSourcePoissonSubtype,CMISSProblemNonlinearSourcePoissonSubtype

  PUBLIC CMISSProblemNoSourceHelmholtzSubtype

  PUBLIC CMISSProblemNoSourceDiffusionSubtype

  PUBLIC CMISSProblemControlSimpleType,CMISSProblemControlFixedLoopType,CMISSProblemControlTimeLoopType, &
    & CMISSProblemControlWhileLoopType

!!
!!==================================================================================================================================
!!
  
CONTAINS

!!==================================================================================================================================

  !>Finalises CMISS.
  SUBROUTINE CMISSFinalise(Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL CMISS_FINALISE(Err,ERROR,*999)

    RETURN
999 CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFinalise

  !
  !================================================================================================================================
  !
  
  !>Initialises CMISS returning a user number to the world coordinate system and region.
  SUBROUTINE CMISSInitialiseNumber(WorldCoordinateUserNumber,WorldRegionUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: WorldCoordinateUserNumber !<On return, the world coordinate system user number.
    INTEGER(INTG), INTENT(OUT) :: WorldRegionUserNumber !<On return, the world region user number.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: WORLD_COORDINATE_SYSTEM
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION

    NULLIFY(WORLD_COORDINATE_SYSTEM)
    NULLIFY(WORLD_REGION)
    CALL CMISS_Initialise(WORLD_REGION,Err,ERROR,*999)
    !CALL CMISS_Initialise(WORLD_COORDINATE_SYSTEM,WORLD_REGION,Err,ERROR,*999)
    WorldCoordinateUserNumber=0
    !WorldCoordinateUserNumber=WORLD_COORDINATE_SYSTEM%USER_NUMBER
    WorldRegionUserNumber=WORLD_REGION%USER_NUMBER

    RETURN
999 CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSInitialiseNumber

  !
  !================================================================================================================================
  !
  
  !>Initialises CMISS returning a pointer to the world coordinate system and region.
  SUBROUTINE CMISSInitialiseObj(WorldCoordinateSystem,WorldRegion,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(OUT) :: WorldCoordinateSystem !<On return, the world coordinate system.
    TYPE(CMISSRegionType), INTENT(OUT) :: WorldRegion !<On return, the world region.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL CMISSCoordinateSystemInitialise(WorldCoordinateSystem)
    CALL CMISSRegionInitialise(WorldRegion)
    CALL CMISS_Initialise(WorldRegion%REGION,Err,ERROR,*999)
    !CALL CMISS_Initialise(WorldCoordinateSystem%COORDINATE_SYSTEM,WorldRegion%REGION,Err,ERROR,*999)

    RETURN
999 CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSInitialiseObj

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSBasisType object.
  SUBROUTINE CMISSBasisTypeFinalise(CMISSBasis,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(OUT) :: CMISSBasis !<The CMISSBasisType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSBasis%BASIS)) CALL BASIS_DESTROY(CMISSBasis%BASIS,Err,ERROR,*999)

    CALL EXITS("CMISSBasisTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSBasisTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSBasisTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSBasisType object.
  SUBROUTINE CMISSBasisTypeInitialise(CMISSBasis,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(OUT) :: CMISSBasis !<The CMISSBasisType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSBasisTypeInitialise",Err,ERROR,*999)

    NULLIFY(CMISSBasis%BASIS)

    CALL EXITS("CMISSBasisTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSBasisTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSBasisTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSBoundaryConditionsType object.
  SUBROUTINE CMISSBoundaryConditionsTypeFinalise(CMISSBoundaryConditions,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(OUT) :: CMISSBoundaryConditions !<The CMISSBoundaryConditionsType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSBoundaryConditionsTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSBoundaryConditions%BOUNDARY_CONDITIONS))  &
      & CALL BOUNDARY_CONDITIONS_DESTROY(CMISSBoundaryConditions%BOUNDARY_CONDITIONS,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsTypeFinalise
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSBoundaryConditionsType object.
  SUBROUTINE CMISSBoundaryConditionsTypeInitialise(CMISSBoundaryConditions,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(OUT) :: CMISSBoundaryConditions !<The CMISSBoundaryConditionsType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSBoundaryConditions%BOUNDARY_CONDITIONS)

    CALL EXITS("CMISSBoundaryConditionsTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSControlLoopType object.
  SUBROUTINE CMISSControlLoopTypeFinalise(CMISSControlLoop,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(OUT) :: CMISSControlLoop !<The CMISSControlLoopType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSControlLoopTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSControlLoop%CONTROL_LOOP))  &
      & CALL CONTROL_LOOP_DESTROY(CMISSControlLoop%CONTROL_LOOP,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeFinalise
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSControlLoopType object.
  SUBROUTINE CMISSControlLoopTypeInitialise(CMISSControlLoop,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(OUT) :: CMISSControlLoop !<The CMISSControlLoopType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSControlLoop%CONTROL_LOOP)

    CALL EXITS("CMISSControlLoopTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSCoordinateSystemType object.
  SUBROUTINE CMISSCoordinateSystemTypeFinalise(CMISSCoordinateSystem,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(OUT) :: CMISSCoordinateSystem !<The CMISSCoordinateSystemType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSCoordinateSystemTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSCoordinateSystem%COORDINATE_SYSTEM))  &
      & CALL COORDINATE_SYSTEM_DESTROY(CMISSCoordinateSystem%COORDINATE_SYSTEM,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeFinalise
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSCoordinateSystemType object.
  SUBROUTINE CMISSCoordinateSystemTypeInitialise(CMISSCoordinateSystem,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(OUT) :: CMISSCoordinateSystem !<The CMISSCoordinateSystemType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSCoordinateSystem%COORDINATE_SYSTEM)

    CALL EXITS("CMISSCoordinateSystemTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSDecompositionType object.
  SUBROUTINE CMISSDecompositionTypeFinalise(CMISSDecomposition,Err)
  
    !Argument variables
    TYPE(CMISSDecompositionType), INTENT(OUT) :: CMISSDecomposition !<The CMISSDecompositionType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSDecompositionTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSDecomposition%DECOMPOSITION))  &
      & CALL DECOMPOSITION_DESTROY(CMISSDecomposition%DECOMPOSITION,Err,ERROR,*999)

    CALL EXITS("CMISSDecompositionTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSDecompositionTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSDecompositionTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSDecompositionTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSDecompositionType object.
  SUBROUTINE CMISSDecompositionTypeInitialise(CMISSDecomposition,Err)
  
    !Argument variables
    TYPE(CMISSDecompositionType), INTENT(OUT) :: CMISSDecomposition !<The CMISSDecompositionType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSDecompositionTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSDecomposition%DECOMPOSITION)

    CALL EXITS("CMISSDecompositionTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSDecompositionTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSDecompositionTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSDecompositionTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSEquationsType object.
  SUBROUTINE CMISSEquationsTypeFinalise(CMISSEquations,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(OUT) :: CMISSEquations !<The CMISSEquationsType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSEquationsTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSEquations%EQUATIONS))  &
      & CALL EQUATIONS_DESTROY(CMISSEquations%EQUATIONS,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSEquationsTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSEquationsTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSEquationsType object.
  SUBROUTINE CMISSEquationsTypeInitialise(CMISSEquations,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(OUT) :: CMISSEquations !<The CMISSEquationsType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSEquationsTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSEquations%EQUATIONS)

    CALL EXITS("CMISSEquationsTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSEquationsTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSEquationsTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSEquationsSetType object.
  SUBROUTINE CMISSEquationsSetTypeFinalise(CMISSEquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: CMISSEquationsSet !<The CMISSEquationsSetType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSEquationsSetTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSEquationsSet%EQUATIONS_SET))  &
      & CALL EQUATIONS_SET_DESTROY(CMISSEquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSEquationsSetTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSEquationsSetTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSEquationsSetType object.
  SUBROUTINE CMISSEquationsSetTypeInitialise(CMISSEquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: CMISSEquationsSet !<The CMISSEquationsSetType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSEquationsSetTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSEquationsSet%EQUATIONS_SET)

    CALL EXITS("CMISSEquationsSetTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSEquationsSetTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSEquationsSetTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSFieldType object.
  SUBROUTINE CMISSFieldTypeFinalise(CMISSField,Err)
  
    !Argument variables
    TYPE(CMISSFieldType), INTENT(OUT) :: CMISSField !<The CMISSFieldType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSFieldTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSField%FIELD))  &
      & CALL FIELD_DESTROY(CMISSField%FIELD,Err,ERROR,*999)

    CALL EXITS("CMISSFieldTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSFieldTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSFieldTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFieldTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSFieldType object.
  SUBROUTINE CMISSFieldTypeInitialise(CMISSField,Err)
  
    !Argument variables
    TYPE(CMISSFieldType), INTENT(OUT) :: CMISSField !<The CMISSFieldType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSFieldTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSField%FIELD)

    CALL EXITS("CMISSFieldTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSFieldTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSFieldTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFieldTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSFieldsType object.
  SUBROUTINE CMISSFieldsTypeFinalise(CMISSFields,Err)
  
    !Argument variables
    TYPE(CMISSFieldsType), INTENT(OUT) :: CMISSFields !<The CMISSFieldsType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSFieldsTypeFinalise",Err,ERROR,*999)
    
    NULLIFY(CMISSFields%FIELDS)

    CALL EXITS("CMISSFieldsTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSFieldsTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSFieldsTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFieldsTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSFieldsType object.
  SUBROUTINE CMISSFieldsTypeInitialise(CMISSFields,Err)
  
    !Argument variables
    TYPE(CMISSFieldsType), INTENT(OUT) :: CMISSFields !<The CMISSFieldsType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSFieldsTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSFields%FIELDS)

    CALL EXITS("CMISSFieldsTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSFieldsTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSFieldsTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFieldsTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSGeneratedMeshType object.
  SUBROUTINE CMISSGeneratedMeshTypeFinalise(CMISSGeneratedMesh,Err)
  
    !Argument variables
    TYPE(CMISSGeneratedMeshType), INTENT(OUT) :: CMISSGeneratedMesh !<The CMISSGeneratedMeshType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSGeneratedMeshTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSGeneratedMesh%GENERATED_MESH))  &
      & CALL GENERATED_MESH_DESTROY(CMISSGeneratedMesh%GENERATED_MESH,Err,ERROR,*999)

    CALL EXITS("CMISSGeneratedMeshTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSGeneratedMeshTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSGeneratedMeshTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSGeneratedMeshTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSGeneratedMeshType object.
  SUBROUTINE CMISSGeneratedMeshTypeInitialise(CMISSGeneratedMesh,Err)
  
    !Argument variables
    TYPE(CMISSGeneratedMeshType), INTENT(OUT) :: CMISSGeneratedMesh !<The CMISSGeneratedMeshType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSGeneratedMeshTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSGeneratedMesh%GENERATED_MESH)

    CALL EXITS("CMISSGeneratedMeshTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSGeneratedMeshTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSGeneratedMeshTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSGeneratedMeshTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSHistoryType object.
  SUBROUTINE CMISSHistoryTypeFinalise(CMISSHistory,Err)
  
    !Argument variables
    TYPE(CMISSHistoryType), INTENT(OUT) :: CMISSHistory !<The CMISSHistoryType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSHistoryTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSHistory%HISTORY))  &
      & CALL HISTORY_DESTROY(CMISSHistory%HISTORY,Err,ERROR,*999)

    CALL EXITS("CMISSHistoryTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSHistoryTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSHistoryTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSHistoryTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSHistoryType object.
  SUBROUTINE CMISSHistoryTypeInitialise(CMISSHistory,Err)
  
    !Argument variables
    TYPE(CMISSHistoryType), INTENT(OUT) :: CMISSHistory !<The CMISSHistoryType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSHistoryTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSHistory%HISTORY)

    CALL EXITS("CMISSHistoryTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSHistoryTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSHistoryTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSHistoryTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSMeshType object.
  SUBROUTINE CMISSMeshTypeFinalise(CMISSMesh,Err)
  
    !Argument variables
    TYPE(CMISSMeshType), INTENT(OUT) :: CMISSMesh !<The CMISSMeshType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSMeshTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSMesh%MESH))  &
      & CALL MESH_DESTROY(CMISSMesh%MESH,Err,ERROR,*999)

    CALL EXITS("CMISSMeshTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSMeshTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSMeshTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSMeshTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSMeshType object.
  SUBROUTINE CMISSMeshTypeInitialise(CMISSMesh,Err)
  
    !Argument variables
    TYPE(CMISSMeshType), INTENT(OUT) :: CMISSMesh !<The CMISSMeshType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSMeshTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSMesh%MESH)

    CALL EXITS("CMISSMeshTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSMeshTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSMeshTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSMeshTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSNodesType object.
  SUBROUTINE CMISSNodesTypeFinalise(CMISSNodes,Err)
  
    !Argument variables
    TYPE(CMISSNodesType), INTENT(OUT) :: CMISSNodes !<The CMISSNodesType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSNodesTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSNodes%NODES))  &
      & CALL NODES_DESTROY(CMISSNodes%NODES,Err,ERROR,*999)

    CALL EXITS("CMISSNodesTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSNodesTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSNodesTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSNodesTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSNodesType object.
  SUBROUTINE CMISSNodesTypeInitialise(CMISSNodes,Err)
  
    !Argument variables
    TYPE(CMISSNodesType), INTENT(OUT) :: CMISSNodes !<The CMISSNodesType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSNodesTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSNodes%NODES)

    CALL EXITS("CMISSNodesTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSNodesTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSNodesTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSNodesTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSProblemType object.
  SUBROUTINE CMISSProblemTypeFinalise(CMISSProblem,Err)
  
    !Argument variables
    TYPE(CMISSProblemType), INTENT(OUT) :: CMISSProblem !<The CMISSProblemType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSProblemTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSProblem%PROBLEM))  &
      & CALL PROBLEM_DESTROY(CMISSProblem%PROBLEM,Err,ERROR,*999)

    CALL EXITS("CMISSProblemTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSProblemTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSProblemTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSProblemTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSProblemType object.
  SUBROUTINE CMISSProblemTypeInitialise(CMISSProblem,Err)
  
    !Argument variables
    TYPE(CMISSProblemType), INTENT(OUT) :: CMISSProblem !<The CMISSProblemType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSProblemTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSProblem%PROBLEM)

    CALL EXITS("CMISSProblemTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSProblemTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSProblemTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSProblemTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSQuadratureType object.
  SUBROUTINE CMISSQuadratureTypeFinalise(CMISSQuadrature,Err)
  
    !Argument variables
    TYPE(CMISSQuadratureType), INTENT(OUT) :: CMISSQuadrature !<The CMISSQuadratureType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSQuadratureTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSQuadrature%QUADRATURE))  &
      & CALL BASIS_QUADRATURE_DESTROY(CMISSQuadrature%QUADRATURE,Err,ERROR,*999)

    CALL EXITS("CMISSQuadratureTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSQuadratureTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSQuadratureTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSQuadratureTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSQuadratureType object.
  SUBROUTINE CMISSQuadratureTypeInitialise(CMISSQuadrature,Err)
  
    !Argument variables
    TYPE(CMISSQuadratureType), INTENT(OUT) :: CMISSQuadrature !<The CMISSQuadratureType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSQuadratureTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSQuadrature%QUADRATURE)

    CALL EXITS("CMISSQuadratureTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSQuadratureTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSQuadratureTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSQuadratureTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSRegionType object.
  SUBROUTINE CMISSRegionTypeFinalise(CMISSRegion,Err)
  
    !Argument variables
    TYPE(CMISSRegionType), INTENT(OUT) :: CMISSRegion !<The CMISSRegionType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSRegionTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSRegion%REGION))  &
      & CALL REGION_DESTROY(CMISSRegion%REGION,Err,ERROR,*999)

    CALL EXITS("CMISSRegionTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSRegionTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSRegionTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSRegionTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSRegionType object.
  SUBROUTINE CMISSRegionTypeInitialise(CMISSRegion,Err)
  
    !Argument variables
    TYPE(CMISSRegionType), INTENT(OUT) :: CMISSRegion !<The CMISSRegionType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSRegionTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSRegion%REGION)

    CALL EXITS("CMISSRegionTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSRegionTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSRegionTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSRegionTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSSolverType object.
  SUBROUTINE CMISSSolverTypeFinalise(CMISSSolver,Err)
  
    !Argument variables
    TYPE(CMISSSolverType), INTENT(OUT) :: CMISSSolver !<The CMISSSolverType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSSolverTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSSolver%SOLVER))  &
      & CALL SOLVER_DESTROY(CMISSSolver%SOLVER,Err,ERROR,*999)

    CALL EXITS("CMISSSolverTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSSolverTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSSolverTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSSolverTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSSolverType object.
  SUBROUTINE CMISSSolverTypeInitialise(CMISSSolver,Err)
  
    !Argument variables
    TYPE(CMISSSolverType), INTENT(OUT) :: CMISSSolver !<The CMISSSolverType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSSolverTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSSolver%SOLVER)

    CALL EXITS("CMISSSolverTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSSolverTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSSolverTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSSolverTypeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSSolverEquationsType object.
  SUBROUTINE CMISSSolverEquationsTypeFinalise(CMISSSolverEquations,Err)
  
    !Argument variables
    TYPE(CMISSSolverEquationsType), INTENT(OUT) :: CMISSSolverEquations !<The CMISSSolverEquationsType object to finalise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSSolverEquationsTypeFinalise",Err,ERROR,*999)
    
    IF(ASSOCIATED(CMISSSolverEquations%SOLVER_EQUATIONS))  &
      & CALL SOLVER_DESTROY(CMISSSolverEquations%SOLVER_EQUATIONS,Err,ERROR,*999)

    CALL EXITS("CMISSSolverEquationsTypeFinalise")
    RETURN
999 CALL ERRORS("CMISSSolverEquationsTypeFinalise",Err,ERROR)
    CALL EXITS("CMISSSolverEquationsTypeFinalise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSSolverEquationsTypeFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSSolverEquationsType object.
  SUBROUTINE CMISSSolverEquationsTypeInitialise(CMISSSolverEquations,Err)
  
    !Argument variables
    TYPE(CMISSSolverEquationsType), INTENT(OUT) :: CMISSSolverEquations !<The CMISSSolverEquationsType object to initialise.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSSolverEquationsTypeInitialise",Err,ERROR,*999)
    
    NULLIFY(CMISSSolverEquations%SOLVER_EQUATIONS)

    CALL EXITS("CMISSSolverEquationsTypeInitialise")
    RETURN
999 CALL ERRORS("CMISSSolverEquationsTypeInitialise",Err,ERROR)
    CALL EXITS("CMISSSolverEquationsTypeInitialise")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSSolverEquationsTypeInitialise

!!==================================================================================================================================
!!
!! ANALYTIC_ANALYSIS_ROUTINES
!!
!!==================================================================================================================================

  !>Output the analytic error analysis for a field specified by a user number compared to the analytic values parameter set.
  SUBROUTINE CMISSAnalyticAnalysisOutputNumber(RegionUserNumber,FieldUserNumber,FileName,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field for analytic error analysis.
    INTEGER(INTG), INTENT(IN) :: FieldUserNumber !<The user number of the field to calculate the analytic error analysis for.
    CHARACTER(LEN=*) :: FileName !<If not empty, the filename to output the analytic analysis to. If empty, the analysis will be output to the standard output.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSAnalyticAnalysisOutputNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL FIELD_USER_NUMBER_FIND(FieldUserNumber,REGION,FIELD,Err,ERROR,*999)
      IF(ASSOCIATED(FIELD)) THEN
        CALL ANALYTIC_ANALYSIS_OUTPUT(FIELD,FileName,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An field with an user number of "//TRIM(NUMBER_TO_VSTRING(FieldUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
    
    CALL EXITS("CMISSAnalyticAnalysisOutputNumber")
    RETURN
999 CALL ERRORS("CMISSAnalyticAnalysisOutputNumber",Err,ERROR)
    CALL EXITS("CMISSAnalyticAnalysisOutputNumber")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSAnalyticAnalysisOutputNumber

  !
  !================================================================================================================================
  !  

  !>Output the analytic error analysis for a field identified by an object compared to the analytic values parameter set.
  SUBROUTINE CMISSAnalyticAnalysisOutputObj(Field,FileName,Err)
  
    !Argument variables
    TYPE(CMISSFieldType), INTENT(IN) :: Field !<The dependent field to calculate the analytic error analysis for.
    CHARACTER(LEN=*) :: FileName !<If not empty, the filename to output the analytic analysis to. If empty, the analysis will be output to the standard output.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSAnalyticAnalysisOutputObj",Err,ERROR,*999)
    
    CALL ANALYTIC_ANALYSIS_OUTPUT(Field%FIELD,FileName,Err,ERROR,*999)

    CALL EXITS("CMISSAnalyticAnalysisOutputObj")
    RETURN
999 CALL ERRORS("CMISSAnalyticAnalysisOutputObj",Err,ERROR)
    CALL EXITS("CMISSAnalyticAnalysisOutputObj")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSAnalyticAnalysisOutputObj

!!==================================================================================================================================
!!
!! BASE_ROUTINES
!!
!!==================================================================================================================================

  !>Sets diagnostics off. \see OPENCMISS::CMISSDiagnosticsSetOn
  SUBROUTINE CMISSDiagnosticsSetOff(Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSDiagnosticsSetOff",Err,ERROR,*999)

    CALL DIAGNOSTICS_SET_OFF(Err,ERROR,*999)

    CALL EXITS("CMISSDiagnosticsSetOff")
    RETURN
999 CALL ERRORS("CMISSDiagnosticsSetOff",Err,ERROR)
    CALL EXITS("CMISSDiagnosticsSetOff")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSDiagnosticsSetOff
  
  !
  !================================================================================================================================
  !
  
  !>Sets diagnostics on \see OPENCMISS::CMISSDiagnosticsSetOff
  SUBROUTINE CMISSDiagnosticsSetOn(DiagType,LevelList,DiagFilename,RoutineList,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: DiagType !<The type of diagnostics to set on \see OPENCMISS_DiagnosticTypes.
    INTEGER(INTG), INTENT(IN) :: LevelList(:) !<The list of diagnostic levels to set on.
    CHARACTER(LEN=*), INTENT(IN) :: DiagFilename !<If present the name of the file to output diagnostic information to. If omitted the diagnostic output is sent to the screen.
    CHARACTER(LEN=*), INTENT(IN) :: RoutineList(:) !<The list of routines to set diagnostics on in.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables
    
    CALL ENTERS("CMISSDiagnosticsSetOn",Err,ERROR,*999)

    CALL DIAGNOSTICS_SET_ON(DiagType,LevelList,DiagFilename,RoutineList,Err,ERROR,*999)

    CALL EXITS("CMISSDiagnosticsSetOn")
    RETURN
999 CALL ERRORS("CMISSDiagnosticsSetOn",Err,ERROR)
    CALL EXITS("CMISSDiagnosticsSetOn")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSDiagnosticsSetOn

  !
  !================================================================================================================================
  !
  
  !>Sets output off \see OPENCMISS::CMISSOutputSetOff
  SUBROUTINE CMISSOutputSetOff(Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSOutputSetOff",Err,ERROR,*999)
    
    CALL OUTPUT_SET_OFF(Err,ERROR,*999)

    CALL EXITS("CMISSOutputSetOff")
    RETURN
999 CALL ERRORS("CMISSOutputSetOff",Err,ERROR)
    CALL EXITS("CMISSOutputSetOff")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSOutputSetOff
  
  !
  !================================================================================================================================
  !
  
  !>Sets output on \see OPENCMISS::CMISSOutputSetOff
  SUBROUTINE CMISSOutputSetOn(EchoFilename,Err)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: EchoFilename !<The filename of the file to echo output to
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSOutputSetOn",Err,ERROR,*999)
    
    CALL OUTPUT_SET_ON(EchoFilename,Err,ERROR,*999)

    CALL EXITS("CMISSOutputSetOn")
    RETURN
999 CALL ERRORS("CMISSOutputSetOn",Err,ERROR)
    CALL EXITS("CMISSOutputSetOn")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSOutputSetOn

  !
  !================================================================================================================================
  !
  
  !>Sets timing off \see OPENCMISS::CMISSTimingSetOn
  SUBROUTINE CMISSTimingSetOff(Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSTimingSetOff",ERR,ERROR,*999)
    
    CALL TIMING_SET_OFF(Err,ERROR,*999)

    CALL EXITS("CMISSTimingSetOff")
    RETURN
999 CALL ERRORS("CMISSTimingSetOff",Err,ERROR)
    CALL EXITS("CMISSTimingSetOff")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSTimingSetOff
  
  !
  !================================================================================================================================
  !
  
  !>Sets timing on \see OPENCMISS::CMISSTimingSetOff
  SUBROUTINE CMISSTimingSetOn(TimingType,TimingSummaryFlag,TimingFilename,RoutineList,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TimingType !<The type of timing to set on \see OPENCMISS_TimingTypes.
    LOGICAL, INTENT(IN) :: TimingSummaryFlag !<.TRUE. if the timing information will be output with subsequent OPENCMISS::CMISSTimingSummaryOutput calls, .FALSE. if the timing information will be output every time the routine exits.
    CHARACTER(LEN=*), INTENT(IN) :: TimingFilename !<If present the name of the file to output timing information to. If omitted the timing output is sent to the screen.
    CHARACTER(LEN=*), INTENT(IN) :: RoutineList(:) !<The list of routines to set timing on in.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSTimingSetOn",Err,ERROR,*999)
    
    CALL TIMING_SET_ON(TimingType,TimingSummaryFlag,TimingFilename,RoutineList,Err,ERROR,*999)

    CALL EXITS("CMISSTimingSetOn")
    RETURN    
999 CALL ERRORS("CMISSTimingSetOn",Err,ERROR)
    CALL EXITS("CMISSTimingSetOn")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSTimingSetOn

  !
  !================================================================================================================================
  !
  
  !>Outputs the timing summary.
  SUBROUTINE CMISSTimingSummaryOutput(Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code     
    !Local variables

    CALL ENTERS("CMISSTimingSummaryOutput",Err,ERROR,*999)
    
    CALL TIMING_SUMMARY_OUTPUT(Err,ERROR,*999)

    CALL EXITS("CMISSTimingSummaryOutput")
    RETURN
999 CALL ERRORS("CMISSTimingSummaryOutput",Err,ERROR)
    CALL EXITS("CMISSTimingSummaryOutput")
    RETURN
    
  END SUBROUTINE CMISSTimingSummaryOutput

!!==================================================================================================================================
!!
!! BASIS_ROUTINES
!!
!!==================================================================================================================================

  !>Returns the collapsed Xi flags of a basis identified by a user number.
  SUBROUTINE CMISSBasisCollapsedXiGetNumber(UserNumber,CollapsedXi,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the collapsed Xi flags for.
    INTEGER(INTG), INTENT(OUT) :: CollapsedXi(:) !<CollapsedXi(ni). On return, the collapsed Xi parameter for the ni'th Xi direction. \see OPENCMISS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSBasisCollapsedXiGetNumber",Err,ERROR,*999)

    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_COLLAPSED_XI_GET(BASIS,CollapsedXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisCollapsedXiGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisCollapsedXiGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisCollapsedXiGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCollapsedXiGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the collapsed Xi flags of a basis identified by an object.
  SUBROUTINE CMISSBasisCollapsedXiGetObj(Basis,CollapsedXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the collapsed Xi flags for.
    INTEGER(INTG), INTENT(OUT) :: CollapsedXi(:) !<CollapsedXi(ni). On return, the collapsed Xi parameter for the ni'th Xi direction. \see OPENCMISS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisCollapsedXiGetObj",Err,ERROR,*999)
    
    CALL BASIS_COLLAPSED_XI_GET(Basis%BASIS,CollapsedXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCollapsedXiGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisCollapsedXiGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisCollapsedXiGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCollapsedXiGetObj

  !
  !================================================================================================================================
  !
 
  !>Sets/changes the collapsed Xi flags of a basis identified by a user number.
  SUBROUTINE CMISSBasisCollapsedXiSetNumber(UserNumber,CollapsedXi,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to set the collapsed Xi flags for.
    INTEGER(INTG), INTENT(IN) :: CollapsedXi(:) !<CollapsedXi(ni). The collapsed Xi parameter for the ni'th Xi direction to set. \see OPENCMISS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisCollapsedXiSetNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_COLLAPSED_XI_SET(BASIS,CollapsedXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisCollapsedXiSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisCollapsedXiSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisCollapsedXiSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCollapsedXiSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the collapsed Xi flags of a basis identified by an object.
  SUBROUTINE CMISSBasisCollapsedXiSetObj(Basis,CollapsedXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to set the collapsed Xi flags for.
    INTEGER(INTG), INTENT(IN) :: CollapsedXi(:) !<CollapsedXi(ni). The collapsed Xi parameter for the ni'th Xi direction to set. \see OPENCMISS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisCollapsedXiSetObj",Err,ERROR,*999)

    CALL BASIS_COLLAPSED_XI_SET(Basis%BASIS,CollapsedXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCollapsedXiSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisCollapsedXiSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisCollapsedXiSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCollapsedXiSetObj

  !
  !================================================================================================================================
  !
  
   !>Finishes the creation of a new basis identified by a user number.
  SUBROUTINE CMISSBasisCreateFinishNumber(UserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisCreateFinishNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,Err,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_CREATE_FINISH(Basis,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF      

    CALL EXITS("CMISSBasisCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSBasisCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSBasisCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCreateFinishNumber

  !
  !================================================================================================================================
  !
  
  !>Finishes the creation of a new basis identified by an object.
  SUBROUTINE CMISSBasisCreateFinishObj(Basis,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSBasisCreateFinishObj",Err,ERROR,*999)

    CALL BASIS_CREATE_FINISH(Basis%BASIS,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSBasisCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSBasisCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCreateFinishObj
  
  !
  !================================================================================================================================
  !
  
  !>Starts the creation of a new basis for a basis identified by a user number.
  SUBROUTINE CMISSBasisCreateStartNumber(UserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to start the creation of.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS

    CALL ENTERS("CMISSBasisCreateStartNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START(UserNumber,BASIS,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSBasisCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSBasisCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCreateStartNumber

  !
  !================================================================================================================================
  !  
 
  !>Starts the creation of a new basis for a basis identified by an object.
  SUBROUTINE CMISSBasisCreateStartObj(UserNumber,Basis,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to start the creation of.
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<On exit, the newly created basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisCreateStartObj",Err,ERROR,*999)
    
    CALL BASIS_CREATE_START(UserNumber,Basis%BASIS,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSBasisCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSBasisCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCreateStartObj
  
  !
  !================================================================================================================================
  !
  
  !>Destroys a basis identified by its basis user number.
  SUBROUTINE CMISSBasisDestroyNumber(UserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisDestroyNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_DESTROY(UserNumber,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
      
    CALL EXITS("CMISSBasisDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSBasisDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSBasisDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisDestroyNumber
  
  !
  !================================================================================================================================
  !
  
  !>Destroys a basis identified by an object.
  SUBROUTINE CMISSBasisDestroyObj(Basis,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisDestroyObj",Err,ERROR,*999)
    
    CALL BASIS_DESTROY(Basis%BASIS,Err,ERROR,*999)

    CALL EXITS("CMISSBasisDestroyObj")
    RETURN
999 CALL ERRORS("CMISSBasisDestroyObj",Err,ERROR)
    CALL EXITS("CMISSBasisDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisDestroyObj
  
  !
  !================================================================================================================================
  !
  
  !>Get the interpolation type in each xi directions for a basis identified by a user number.
  SUBROUTINE CMISSBasisInterpolationXiGetNumber(UserNumber,InterpolationXi,Err)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the interpolation xi for.
    INTEGER(INTG), INTENT(OUT) :: InterpolationXi(:) !<On return, the interpolation xi parameters for each Xi direction \see OPENCMISS_InterpolationSpecifications.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("CMISSBasisInterpolationXiGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_INTERPOLATION_XI_GET(BASIS,InterpolationXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisInterpolationXiGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisInterpolationXiGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisInterpolationXiGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisInterpolationXiGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Get the interpolation type in each xi directions for a basis indentified by an object. 
  SUBROUTINE CMISSBasisInterpolationXiGetObj(Basis,InterpolationXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the interpolation xi for.
    INTEGER(INTG), INTENT(OUT) :: InterpolationXi(:) !<On return, the interpolation xi parameters for each Xi direction \see OPENCMISS_InterpolationSpecifications.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisInterpolationXiGetObj",Err,ERROR,*999)
    
    CALL BASIS_INTERPOLATION_XI_GET(Basis%BASIS,InterpolationXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisInterpolationXiGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisInterpolationXiGetObj",Err,ERROR)
    CALL CMISS_HANDL_ERROR(Err,ERROR)
    CALL EXITS("CMISSBasisInterpolationXiGetObj")
    RETURN
    
  END SUBROUTINE CMISSBasisInterpolationXiGetObj
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the interpolation type in each xi directions for a basis identified by a user number.
  SUBROUTINE CMISSBasisInterpolationXiSetNumber(UserNumber,InterpolationXi,Err)
  
    !Argument variables
     INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the interpolation xi for.
     INTEGER(INTG), INTENT(IN) :: InterpolationXi(:) !<The interpolation xi parameters for each Xi direction \see OPENCMISS_InterpolationSpecifications.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisInterpolationXiSetNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_INTERPOLATION_XI_SET(BASIS,InterpolationXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF      

    CALL EXITS("CMISSBasisInterpolationXiSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisInterpolationXiSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisInterpolationXiSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisInterpolationXiSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the interpolation type in each xi directions for a basis indentified by an object.
  SUBROUTINE CMISSBasisInterpolationXiSetObj(Basis,InterpolationXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the interpolation xi for.
    INTEGER(INTG), INTENT(IN) :: InterpolationXi(:) !<The interpolation xi parameters for each Xi direction \see OPENCMISS_InterpolationSpecifications.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisInterpolationXiSetObj",Err,ERROR,*999)
    
    CALL BASIS_INTERPOLATION_XI_SET(Basis%BASIS,InterpolationXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisInterpolationXiSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisInterpolationXiSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisInterpolationXiSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisInterpolationXiSetObj
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of local nodes in a basis identified by a user number.
  SUBROUTINE CMISSBasisNumberOfLocalNodesGetNumber(UserNumber,NumberOfLocalNodes,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the interpolation xi for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfLocalNodes !<On return, the number of local nodes in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisNumberOfLocalNodesGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_NUMBER_OF_LOCAL_NODES_GET(BASIS,NumberOfLocalNodes,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisNumberOfLocalNodesGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfLocalNodesGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfLocalNodesGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfLocalNodesGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of local nodes in a basis identified by an object.
  SUBROUTINE CMISSBasisNumberOfLocalNodesGetObj(Basis,NumberOfLocalNodes,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the number of local nodes for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfLocalNodes !<On return, the number of local nodes in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisNumberOfLocalNodesGetObj",Err,ERROR,*999)

    CALL BASIS_NUMBER_OF_LOCAL_NODES_GET(Basis%BASIS,NumberOfLocalNodes,Err,ERROR,*999)

    CALL EXITS("CMISSBasisNumberOfLocalNodesGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfLocalNodesGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfLocalNodesGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfLocalNodesGetObj
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Xi directions in a basis identified by a user number.
  SUBROUTINE CMISSBasisNumberOfXiGetNumber(UserNumber,NumberOfXi,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the number xi for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfXi !<On return, the number of xi directions in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisNumberOfXiGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_NUMBER_OF_XI_GET(BASIS,NumberOfXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisNumberOfXiGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfXiGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfXiGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfXiGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Xi directions in a basis identified by an object.
  SUBROUTINE CMISSBasisNumberOfXiGetObj(Basis,NumberOfXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !The basis to get the number of xi directions for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfXi !<On return, the number of xi directions in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisNumberOfXiGetObj",Err,ERROR,*999)
    
    CALL BASIS_NUMBER_OF_XI_GET(Basis%BASIS,NumberOfXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisNumberOfXiGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfXiGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfXiGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfXiGetObj
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the number of Xi directions in a basis identified by a user number.
  SUBROUTINE CMISSBasisNumberOfXiSetNumber(UserNumber,NumberOfXi,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to set the number xi for.
    INTEGER(INTG), INTENT(IN) :: NumberOfXi !<The number of xi directions in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisNumberOfXiSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_NUMBER_OF_XI_SET(BASIS,NumberOfXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisNumberOfXiSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfXiSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfXiSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfXiSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the number of Xi directions in a basis identified by an object.
  SUBROUTINE CMISSBasisNumberOfXiSetObj(Basis,NumberOfXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to set the number of xi directions for.
    INTEGER(INTG), INTENT(IN) :: NumberOfXi !<The number of xi directions in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisNumberOfXiSetObj",Err,ERROR,*999)
    
    CALL BASIS_NUMBER_OF_XI_SET(Basis%BASIS,NumberOfXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisNumberOfXiSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfXiSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfXiSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfXiSetObj
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Gauss points in each Xi directions for a basis quadrature identified by a user number.
  SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetNumber(UserNumber,NumberOfGaussXi,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the number of Gauss Xi for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfGaussXi(:) !<On return, the number of Gauss points in each Xi directions in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET(BASIS,NumberOfGaussXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number Gauss points in each Xi directions for a basis quadrature identified by an object.
  SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetObj(Basis,NumberOfGaussXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the number of Gauss Xi for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfGaussXi(:) !<On return, the number of Gauss points in each Xi directions in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiGetObj",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET(Basis%BASIS,NumberOfGaussXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetObj
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the number of Gauss points in each Xi directions for a basis quadrature identified by a user number.
  SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiSetNumber(UserNumber,NumberOfGaussXi,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to set the number of Gauss Xi for.
    INTEGER(INTG), INTENT(IN) :: NumberOfGaussXi(:) !<The number of Gauss points in each Xi directions in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(BASIS,NumberOfGaussXi,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
    
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets the number Gauss points in each Xi directions for a basis quadrature identified by an object.
  SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiSetObj(Basis,NumberOfGaussXi,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to get the number of Gauss Xi for.
    INTEGER(INTG), INTENT(IN) :: NumberOfGaussXi(:) !<The number of Gauss points in each Xi directions in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiSetObj",Err,ERROR,*999)

    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(Basis%BASIS,NumberOfGaussXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiSetObj
  
  !
  !================================================================================================================================
  !
  
  !>Returns the order of quadrature a basis quadrature identified by a user number.
  SUBROUTINE CMISSBasisQuadratureOrderGetNumber(UserNumber,QuadratureOrder,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the quadrature order for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureOrder !<On return, the order of quadrature in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisQuadratureOrderGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_QUADRATURE_ORDER_GET(BASIS,QuadratureOrder,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisQuadratureOrderGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureOrderGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureOrderGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureOrderGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the the order of quadrature for a basis quadrature identified by an object.
  SUBROUTINE CMISSBasisQuadratureOrderGetObj(Basis,QuadratureOrder,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the quadrature order for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureOrder !<On return, the order of quadrature in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureOrderGetObj",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_ORDER_GET(Basis%BASIS,QuadratureOrder,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureOrderGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureOrderGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureOrderGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureOrderGetObj

  !
  !================================================================================================================================
  !
    
  !>Sets/changes the order of quadrature a basis quadrature identified by a user number.
  SUBROUTINE CMISSBasisQuadratureOrderSetNumber(UserNumber,QuadratureOrder,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to set the quadrature order for.
    INTEGER(INTG), INTENT(IN) :: QuadratureOrder !<The order of quadrature in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisQuadratureOrderSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_QUADRATURE_ORDER_SET(BASIS,QuadratureOrder,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisQuadratureOrderSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureOrderSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureOrderSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureOrderSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the the order of quadrature for a basis quadrature identified by an object.
  SUBROUTINE CMISSBasisQuadratureOrderSetObj(Basis,QuadratureOrder,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to set the quadrature order for.
    INTEGER(INTG), INTENT(IN) :: QuadratureOrder !<The order of quadrature in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureOrderSetObj",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_ORDER_SET(Basis%BASIS,QuadratureOrder,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureOrderSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureOrderSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureOrderSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureOrderSetObj
  
  !
  !================================================================================================================================
  !
  
  !>Returns the type of quadrature a basis quadrature identified by a user number.
  SUBROUTINE CMISSBasisQuadratureTypeGetNumber(UserNumber,QuadratureType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the quadrature type for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureType !<On return, the type of quadrature in the specified basis. \see OPENCMISS_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisQuadratureTypeGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_QUADRATURE_TYPE_GET(BASIS,QuadratureType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisQuadratureTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureTypeGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the the type of quadrature for a basis quadrature identified by an object.
  SUBROUTINE CMISSBasisQuadratureTypeGetObj(Basis,QuadratureType,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the quadrature order for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureType !<On return, the type of quadrature in the specified basis. \see OPENCMISS_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureTypeGetObj",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_TYPE_GET(Basis%BASIS,QuadratureType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureTypeGetObj

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the type of quadrature a basis quadrature identified by a user number.
  SUBROUTINE CMISSBasisQuadratureTypeSetNumber(UserNumber,QuadratureType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the quadrature type for.
    INTEGER(INTG), INTENT(IN) :: QuadratureType !<The type of quadrature in the specified basis to set. \see OPENCMISS_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisQuadratureTypeSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_QUADRATURE_TYPE_SET(BASIS,QuadratureType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisQuadratureTypeSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureTypeSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureTypeSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureTypeSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the the type of quadrature for a basis quadrature identified by an object.
  SUBROUTINE CMISSBasisQuadratureTypeSetObj(Basis,QuadratureType,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to get the quadrature type for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureType !<The type of quadrature in the specified basis to set. \see OPENCMISS_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureTypeSetObj",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_TYPE_SET(Basis%BASIS,QuadratureType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureTypeSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureTypeSetObj

  !
  !================================================================================================================================
  !
  
  !>Returns the type of a basis identified by a user number.
  SUBROUTINE CMISSBasisTypeGetNumber(UserNumber,BasisType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to get the type for.
    INTEGER(INTG), INTENT(OUT) :: BasisType !<On return, the type of the specified basis. \see OPENCMISS_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisTypeGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_TYPE_GET(BASIS,BasisType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the type of a basis identified by an object.
  SUBROUTINE CMISSBasisTypeGetObj(Basis,BasisType,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(IN) :: Basis !<The basis to get the type for.
    INTEGER(INTG), INTENT(OUT) :: BasisType !<On return, the type of the specified basis. \see OPENCMISS_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisTypeGetObj",Err,ERROR,*999)
    
    CALL BASIS_TYPE_GET(Basis%BASIS,BasisType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSBasisTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSBasisTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeGetObj

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the type of a basis identified by a user number.
  SUBROUTINE CMISSBasisTypeSetNumber(UserNumber,BasisType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to set the type for.
    INTEGER(INTG), INTENT(IN) :: BasisType !<The type of the specified basis to set. \see OPENCMISS_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBasisTypeSetNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN
      CALL BASIS_TYPE_SET(BASIS,BasisType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A basis with an user number of "//TRIM(NUMBER_TO_VSTRING(UserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBasisTypeSetNumber")
    RETURN
999 CALL ERRORS("CMISSBasisTypeSetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisTypeSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeSetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the type of a basis identified by an object.
  SUBROUTINE CMISSBasisTypeSetObj(Basis,BasisType,Err)
  
    !Argument variables
    TYPE(CMISSBasisType), INTENT(INOUT) :: Basis !<The basis to set the type for.
    INTEGER(INTG), INTENT(IN) :: BasisType !<The type of the specified basis to set. \see OPENCMISS_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisTypeSetObj",Err,ERROR,*999)
    
    CALL BASIS_TYPE_SET(Basis%BASIS,BasisType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSBasisTypeSetObj",Err,ERROR)
    CALL EXITS("CMISSBasisTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeSetObj


!!==================================================================================================================================
!!
!! BOUNDARY_CONDITIONS_ROUTINES
!!
!!==================================================================================================================================
  
  !>Destroys the boundary conditions for an equations set identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to destroy the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the boundary conditions for.
     INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsDestroyNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_DESTROY(BOUNDARY_CONDITIONS,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsDestroyNumber

  !
  !================================================================================================================================
  !
    
  !>Destroys boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsDestroyObj(BoundaryConditions,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(INOUT) :: BoundaryConditions !<The boundary conditions to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsDestroyObj",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_DESTROY(BoundaryConditions%BOUNDARY_CONDITIONS,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsDestroyObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsDestroyObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsDestroyObj

  !
  !================================================================================================================================
  !
  
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsAddConstantNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,ComponentNumber, &
    &  Condition,Value,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to add the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsAddConstantNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
       IF(ASSOCIATED(EQUATIONS_SET)) THEN
         CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
         CALL BOUNDARY_CONDITIONS_ADD_CONSTANT(BOUNDARY_CONDITIONS,VariableType,ComponentNumber,Condition,Value,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsAddConstantNumber")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddConstantNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddConstantNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddConstantNumber
  
  !
  !================================================================================================================================
  !
  
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant for boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsAddConstantObj(BoundaryConditions,VariableType,ComponentNumber,Condition,Value,Err)
    
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(IN) :: BoundaryConditions !<The boundary conditions to add the constant to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSBoundaryConditionsAddConstantObj",Err,ERROR,*999)
  
    CALL BOUNDARY_CONDITIONS_ADD_CONSTANT(BoundaryConditions%BOUNDARY_CONDITIONS,VariableType,ComponentNumber,Condition,Value, &
      & Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsAddConstantObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddConstantObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddConstantObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddConstantObj
 
  !
  !================================================================================================================================
  !
  

  !>Sets the value of the specified constant as a boundary condition on the specified constant for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsSetConstantNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsSetConstantNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_SET_CONSTANT(BOUNDARY_CONDITIONS,VariableType,ComponentNumber,Condition,Value,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsSetConstantNumber")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetConstantNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetConstantNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetConstantNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets the value of the specified constant and sets this as a boundary condition on the specified constant for boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsSetConstantObj(BoundaryConditions,VariableType,ComponentNumber,Condition,Value,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(IN) :: BoundaryConditions !<The boundary conditions to set the constant to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSBoundaryConditionsSetConstantObj",ERR,ERROR,*999)
  
    CALL BOUNDARY_CONDITIONS_SET_CONSTANT(BoundaryConditions%BOUNDARY_CONDITIONS,VariableType,ComponentNumber,Condition,Value, &
      & Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsSetConstantObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetConstantObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetConstantObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetConstantObj

  !
  !================================================================================================================================
  !
  
  !>Adds the value to the specified element and sets this as a boundary condition on the specified element for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsAddElementNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,ElementUserNumber, &
    & ComponentNumber,Condition,Value,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ElementUserNumber !<The user number of the element to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to add the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsAddElementNumber",Err,ERROR,*999)

    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_ADD_ELEMENT(BOUNDARY_CONDITIONS,VariableType,ElementUserNumber,ComponentNumber,Condition,Value, &
          & Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsAddElementNumber")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddElementNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddElementNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddElementNumber
  
  !
  !================================================================================================================================
  !
  
  !>Adds to the value of the specified element and sets this as a boundary condition on the specified element for boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsAddElementObj(BoundaryConditions,VariableType,ElementUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(IN) :: BoundaryConditions !<The boundary conditions to add the element to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ElementUserNumber !<The user number of the element to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsAddElementObj",Err,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_ADD_ELEMENT(BoundaryConditions%BOUNDARY_CONDITIONS,VariableType,ElementUserNumber,ComponentNumber, &
      & Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsAddElementObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddElementObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddElementObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddElementObj
 
  !
  !================================================================================================================================
  !

  !>Sets the value of the specified element as a boundary condition on the specified element for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsSetElementNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,ElementUserNumber, &
    & ComponentNumber,Condition,Value,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ElementUserNumber !<The user number of the element to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsSetElementNumber",Err,ERROR,*999)

    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_SET_ELEMENT(BOUNDARY_CONDITIONS,VariableType,ElementUserNumber,ComponentNumber,Condition,Value, &
          & Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsSetElementNumber")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetElementNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetElementNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetElementNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets the value of the specified element and sets this as a boundary condition on the specified elements for boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsSetElementObj(BoundaryConditions,VariableType,ElementUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(IN) :: BoundaryConditions !<The boundary conditions to set the element to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ElementUserNumber !<The user number of the element to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsSetElementObj",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_SET_ELEMENT(BoundaryConditions%BOUNDARY_CONDITIONS,VariableType,ElementUserNumber,ComponentNumber, &
      & Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsSetElementObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetElementObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetElementObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetElementObj

  !
  !================================================================================================================================
  !
  
  !>Adds the value to the specified node and sets this as a boundary condition on the specified node for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsAddNodeNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,DerivativeNumber, &
    & NodeUserNumber,ComponentNumber,Condition,Value,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: DerivativeNumber !<The user number of the node derivative to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: NodeUserNumber !<The user number of the node to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to add the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsAddNodeNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_ADD_NODE(BOUNDARY_CONDITIONS,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
          & Condition,Value,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsAddNodeNumber")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddNodeNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddNodeNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddNodeNumber
  
  !
  !================================================================================================================================
  !
  
  !>Adds to the value of the specified node and sets this as a boundary condition on the specified node for boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsAddNodeObj(BoundaryConditions,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(IN) :: BoundaryConditions !<The boundary conditions to add the node to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: DerivativeNumber !<The user number of the node derivative to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: NodeUserNumber !<The user number of the node to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsAddNodeObj",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_ADD_NODE(BoundaryConditions%BOUNDARY_CONDITIONS,VariableType,DerivativeNumber,NodeUserNumber, &
      & ComponentNumber,Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsAddNodeObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddNodeObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddNodeObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddNodeObj
 
  !
  !================================================================================================================================
  !

  !>Sets the value of the specified node as a boundary condition on the specified node for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsSetNodeNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,DerivativeNumber, &
    & NodeUserNumber,ComponentNumber,Condition,Value,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: DerivativeNumber !<The user number of the node derivative to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: NodeUserNumber !<The user number of the node to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSBoundaryConditionsSetNodeNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
          & Condition,Value,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSBoundaryConditionsSetNodeNumber")
    RETURN    
999 CALL ERRORS("CMISSBoundaryConditionsSetNodeNumber",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetNodeNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetNodeNumber
  
  !
  !================================================================================================================================
  !
  
  !>Sets the value of the specified node and sets this as a boundary condition on the specified node for boundary conditions identified by an object.
  SUBROUTINE CMISSBoundaryConditionsSetNodeObj(BoundaryConditions,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(CMISSBoundaryConditionsType), INTENT(IN) :: BoundaryConditions !<The boundary conditions to set the node to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: DerivativeNumber !<The user number of the node derivative to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: NodeUserNumber !<The user number of the node to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSBoundaryConditionsSetNodeObj",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_SET_NODE(BoundaryConditions%BOUNDARY_CONDITIONS,VariableType,DerivativeNumber,NodeUserNumber, &
      & ComponentNumber,Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsSetNodeObj")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetNodeObj",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetNodeObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetNodeObj

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for an equations set identified by a user number. 
  SUBROUTINE CMISSEquationsSetBoundaryConditionsGetNumber(RegionUserNumber,EquationsSetUserNumber,BoundaryConditions,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to get the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the boundary conditions for.
    TYPE(CMISSBoundaryConditionsType), INTENT(OUT) :: BoundaryConditions !<On return, The boundary conditions for the specified equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetBoundaryConditionsGetNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
       IF(ASSOCIATED(EQUATIONS_SET)) THEN
         CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BoundaryConditions%BOUNDARY_CONDITIONS,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetBoundaryConditionsGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Gets the boundary conditions for an equations set identified by a user number. 
  SUBROUTINE CMISSEquationsSetBoundaryConditionsGetObj(EquationsSet,BoundaryConditions,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(IN) :: EquationsSet !<The equations set to get the boundary conditions for.
    TYPE(CMISSBoundaryConditionsType), INTENT(OUT) :: BoundaryConditions !<On return, the boundary conditions for the specified equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsGetObj",ERR,ERROR,*999)

    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EquationsSet%EQUATIONS_SET,BoundaryConditions%BOUNDARY_CONDITIONS,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetBoundaryConditionsGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsGetObj

!!==================================================================================================================================
!!
!! CONTROL_LOOP_ROUTINES
!!
!!==================================================================================================================================

  !>Gets the current time parameters for a time control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopCurrentTimesGetNumber0(ProblemUserNumber,ControlLoopIdentifier,CurrentTime,TimeIncrement,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time of the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the current time increment of the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopCurrentTimesGetNumber0",Err,ERROR,*999)
 
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CurrentTime,TimeIncrement,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopCurrentTimesGetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopCurrentTimesGetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopCurrentTimesGetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopCurrentTimesGetNumber0

  !
  !================================================================================================================================
  !  
  
  !>Gets the current time parameters for a time control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopCurrentTimesGetNumber1(ProblemUserNumber,ControlLoopIdentifiers,CurrentTime,TimeIncrement,Err)
  
    !Argument variables
     INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time of the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the current time increment of the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopCurrentTimesGetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CurrentTime,TimeIncrement,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopCurrentTimesGetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopCurrentTimesGetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopCurrentTimesGetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopCurrentTimesGetNumber1

  !
  !================================================================================================================================
  !  
  
  !>Gets the current time parameters for a time control loop identified by an object.
  SUBROUTINE CMISSControlLoopCurrentTimesGetObj(ControlLoop,CurrentTime,TimeIncrement,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(IN) :: ControlLoop !<The control loop to get the current times for.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time of the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the current time increment of the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopCurrentTimesGetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_CURRENT_TIMES_GET(ControlLoop%CONTROL_LOOP,CurrentTime,TimeIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopCurrentTimesGetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopCurrentTimesGetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopCurrentTimesGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopCurrentTimesGetObj

  !
  !================================================================================================================================
  !  
  
  !>Destroys a control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopDestroyNumber0(ProblemUserNumber,ControlLoopIdentifier,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to destroy the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopDestroyNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_DESTROY(CONTROL_LOOP,Err,ERROR,*999)
   ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
    
    CALL EXITS("CMISSControlLoopDestroyNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopDestroyNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopDestroyNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopDestroyNumber0

  !
  !================================================================================================================================
  !  
  
  !>Destroys a control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopDestroyNumber1(ProblemUserNumber,ControlLoopIdentifiers,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to destroy the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopDestroyNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_DESTROY(CONTROL_LOOP,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopDestroyNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopDestroyNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopDestroyNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopDestroyNumber1
  
  !
  !================================================================================================================================
  !  
  
  !>Destroys a control loop identified by an object.
  SUBROUTINE CMISSControlLoopDestroyObj(ControlLoop,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(INOUT) :: ControlLoop !<The control loop to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopDestroyObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_DESTROY(ControlLoop%CONTROL_LOOP,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopDestroyObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopDestroyObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopDestroyObj
 
  !
  !================================================================================================================================
  !  
  
  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root identified by user numbers.
  SUBROUTINE CMISSControlLoopGetNumber00(ProblemUserNumber,ControlLoopRootIdentifier,ControlLoopIdentifier,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifier !<The root control loop identifier.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    TYPE(CMISSControlLoopType), INTENT(OUT) :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopGetNumber00",Err,ERROR,*999)
    
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifier,ROOT_CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifier,ControlLoop%CONTROL_LOOP,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopGetNumber00")
    RETURN    
999 CALL ERRORS("CMISSControlLoopGetNumber00",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetNumber00")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetNumber00

  !
  !================================================================================================================================
  !  
  
  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root identified by user numbers.
  SUBROUTINE CMISSControlLoopGetNumber10(ProblemUserNumber,ControlLoopRootIdentifiers,ControlLoopIdentifier,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifiers(:) !<The root control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    TYPE(CMISSControlLoopType), INTENT(OUT) :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopGetNumber10",Err,ERROR,*999)
 
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifiers,ROOT_CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifier,ControlLoop%CONTROL_LOOP,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopGetNumber10")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetNumber10",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetNumber10")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetNumber10

  !
  !================================================================================================================================
  !  
  
  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root identified by user numbers.
  SUBROUTINE CMISSControlLoopGetNumber01(ProblemUserNumber,ControlLoopRootIdentifier,ControlLoopIdentifiers,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifier !<The root control loop identifier.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    TYPE(CMISSControlLoopType), INTENT(OUT) :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopGetNumber01",Err,ERROR,*999)
    
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifier,ROOT_CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifiers,ControlLoop%CONTROL_LOOP,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopGetNumber01")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetNumber01",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetNumber01")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetNumber01

  !
  !================================================================================================================================
  !  
  
  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root identified by user numbers.
  SUBROUTINE CMISSControlLoopGetNumber11(ProblemUserNumber,ControlLoopRootIdentifiers,ControlLoopIdentifiers,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifiers(:) !<The root control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    TYPE(CMISSControlLoopType), INTENT(OUT) :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopGetNumber11",Err,ERROR,*999)
    
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifiers,ROOT_CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifiers,ControlLoop%CONTROL_LOOP,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXTIS("CMISSControlLoopGetNumber11")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetNumber11",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetNumber11")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetNumber11

  !
  !================================================================================================================================
  !  
  
  !>Destroys a control loop identified by an object.
  SUBROUTINE CMISSControlLoopGetObj0(ControlLoopRoot,ControlLoopIdentifier,ControlLoop,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(IN) :: ControlLoopRoot !<The root control loop.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    TYPE(CMISSControlLoopType), INTENT(OUT) :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopGetObj0",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_GET(ControlLoopRoot%CONTROL_LOOP,ControlLoopIdentifier,ControlLoop%CONTROL_LOOP,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopGetObj0")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetObj0",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetObj0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetObj0
  
  !
  !================================================================================================================================
  !  
  
  !>Destroys a control loop identified by an object.
  SUBROUTINE CMISSControlLoopGetObj1(ControlLoopRoot,ControlLoopIdentifiers,ControlLoop,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(IN) :: ControlLoopRoot !<The root control loop.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    TYPE(CMISSControlLoopType), INTENT(OUT) :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopGetObj1",Err,ERROR,*999)

    CALL CONTROL_LOOP_GET(ControlLoopRoot%CONTROL_LOOP,ControlLoopIdentifiers,ControlLoop%CONTROL_LOOP,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopGetObj1")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetObj1",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetObj1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetObj1

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the iteration parameters for a fixed control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopIterationsSetNumber0(ProblemUserNumber,ControlLoopIdentifier,StartIteration,StopIteration, &
    & IterationIncrement,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the iteration parameters for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    INTEGER(INTG), INTENT(IN) :: StartIteration !<The start iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: StopIteration !<The stop iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: IterationIncrement !<The iteration increment of the fixed control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopIterationsSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_ITERATIONS_SET(CONTROL_LOOP,StartIteration,StopIteration,IterationIncrement,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopIterationsSetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopIterationsSetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopIterationsSetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopIterationsSetNumber0

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the iteration parameters for a fixed control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopIterationsSetNumber1(ProblemUserNumber,ControlLoopIdentifiers,StartIteration,StopIteration, &
    & IterationIncrement,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the iteration parameters for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: StartIteration !<The start iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: StopIteration !<The stop iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: IterationIncrement !<The iteration increment of the fixed control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopIterationsSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_ITERATIONS_SET(CONTROL_LOOP,StartIteration,StopIteration,IterationIncrement,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopIterationsSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopIterationsSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopIterationsSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopIterationsSetNumber1

  !================================================================================================================================
  !  
  
  !>Sets/changes the iteration parameters for a fixed control loop identified by an object.
  SUBROUTINE CMISSControlLoopIterationsSetObj(ControlLoop,StartIteration,StopIteration,IterationIncrement,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(INOUT) :: ControlLoop !<The control loop to set the iteration parameters for.
    INTEGER(INTG), INTENT(IN) :: StartIteration !<The start iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: StopIteration !<The stop iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: IterationIncrement !<The iteration increment of the fixed control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopIterationsSetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_ITERATIONS_SET(ControlLoop%CONTROL_LOOP,StartIteration,StopIteration,IterationIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopIterationsSetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopIterationsSetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopIterationsSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopIterationsSetObj

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the maximum iterations for a while control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopMaximumIterationsSetNumber0(ProblemUserNumber,ControlLoopIdentifier,MaximumIterations,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the maximum iterations for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    INTEGER(INTG), INTENT(IN) :: MaximumIterations !<The maximum iterations of the while control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopMaximumIterationsSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MaximumIterations,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopMaximumIterationsSetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopMaximumIterationsSetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopMaximumIterationsSetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopMaximumIterationsSetNumber0

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the maximum iterations for a while control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopMaximumIterationsSetNumber1(ProblemUserNumber,ControlLoopIdentifiers,MaximumIterations,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the maximum iterations for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: MaximumIterations !<The maximum iterations of the while control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopMaximumIterationsSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MaximumIterations,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
    
    CALL EXITS("CMISSControlLoopMaximumIterationsSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopMaximumIterationsSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopMaximumIterationsSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopMaximumIterationsSetNumber1

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the maximum iterations for a while control loop identified by an object.
  SUBROUTINE CMISSControlLoopMaximumIterationsSetObj(ControlLoop,MaximumIterations,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(INOUT) :: ControlLoop !<The control loop to set the maximum iterations for.
    INTEGER(INTG), INTENT(IN) :: MaximumIterations !<The maximum iterations of the while control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopMaximumIterationsSetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(ControlLoop%CONTROL_LOOP,MaximumIterations,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopMaximumIterationsSetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopMaximumIterationsSetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopMaximumIterationsSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopMaximumIterationsSetObj

  !
  !================================================================================================================================
  !  
  
  !>Returns the number of sub-control loops for a control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetNumber0(ProblemUserNumber,ControlLoopIdentifier,NumberOfSubLoops,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the number of sub loops for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    INTEGER(INTG), INTENT(OUT) :: NumberOfSubLoops !<On return, the number of sub loops for the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsGetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsGetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetNumber0

  !
  !================================================================================================================================
  !  
  
  !>Returns the number of sub-control loops for a control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetNumber1(ProblemUserNumber,ControlLoopIdentifiers,NumberOfSubLoops,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the number of sub loops for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    INTEGER(INTG), INTENT(OUT) :: NumberOfSubLoops !<On return, the number of sub loops for the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsGetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsGetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetNumber1

  !
  !================================================================================================================================
  !  
  
  !>Returns the number of sub-control loops for a control loop identified by an object.
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetObj(ControlLoop,NumberOfSubLoops,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(IN) :: ControlLoop !<The control loop to get the number of sub loops for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfSubLoops !<On return, the number of sub loops for the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsGetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(ControlLoop%CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsGetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetObj

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the number of sub-control loops for a control loop identified by user numbers. \todo is this really public???
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetNumber0(ProblemUserNumber,ControlLoopIdentifier,NumberOfSubLoops,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the number of sub loops for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    INTEGER(INTG), INTENT(IN) :: NumberOfSubLoops !<The number of sub loops for the specified control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsSetNumber",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsSetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetNumber0

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the number of sub-control loops for a control loop identified by user numbers. \todo is this really public???
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetNumber1(ProblemUserNumber,ControlLoopIdentifiers,NumberOfSubLoops,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the number of sub loops for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: NumberOfSubLoops !<The number of sub loops for the specified control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
    
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetNumber1

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the number of sub-control loops for a control loop identified by an object. \todo is this really public???
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetObj(ControlLoop,NumberOfSubLoops,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(INOUT) :: ControlLoop !<The control loop to set the number of sub loops for.
    INTEGER(INTG), INTENT(IN) :: NumberOfSubLoops !<The number of sub loops for the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsSetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(ControlLoop%CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsSetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetObj

  !
  !================================================================================================================================
  !  
  
  !>Returns the time parameters for a time control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopTimesGetNumber0(ProblemUserNumber,ControlLoopIdentifier,StartTime,StopTime,TimeIncrement, &
    & CurrentTime,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the time parameters for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    REAL(DP), INTENT(OUT) :: StartTime !<On return, the start time for the time control loop.
    REAL(DP), INTENT(OUT) :: StopTime !<On return, the stop time for the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the time increment for the time control loop.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopTimesGetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,CurrentTime,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopTimesGetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesGetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesGetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesGetNumber0

  !
  !================================================================================================================================
  !  
 
  !>Returns the time parameters for a time control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopTimesGetNumber1(ProblemUserNumber,ControlLoopIdentifiers,StartTime,StopTime,TimeIncrement, &
    & CurrentTime,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the time parameters for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifier.
    REAL(DP), INTENT(OUT) :: StartTime !<On return, the start time for the time control loop.
    REAL(DP), INTENT(OUT) :: StopTime !<On return, the stop time for the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the time increment for the time control loop.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopTimesGetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,CurrentTime,Err,ERROR,*999)
   ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF
    
    CALL EXITS("CMISSControlLoopTimesGetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesGetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesGetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesGetNumber1

  !
  !================================================================================================================================
  !  
 
  !>Returns the time parameters for a time control loop identified by an object.
  SUBROUTINE CMISSControlLoopTimesGetObj(ControlLoop,StartTime,StopTime,TimeIncrement,CurrentTime,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(IN) :: ControlLoop !<The control loop to get the times for.
    REAL(DP), INTENT(OUT) :: StartTime !<On return, the start time for the time control loop.
    REAL(DP), INTENT(OUT) :: StopTime !<On return, the stop time for the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the time increment for the time control loop.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTimesGetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_TIMES_GET(ControlLoop%CONTROL_LOOP,StartTime,StopTime,TimeIncrement,CurrentTime,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTimesGetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesGetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesGetObj

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the time parameters for a time control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopTimesSetNumber0(ProblemUserNumber,ControlLoopIdentifier,StartTime,StopTime,TimeIncrement,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the time parameters for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    REAL(DP), INTENT(IN) :: StartTime !<The start time for the time control loop to set.
    REAL(DP), INTENT(IN) :: StopTime !<The stop time for the time control loop to set.
    REAL(DP), INTENT(IN) :: TimeIncrement !<The time increment for the time control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopTimesSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopTimesSetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesSetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesSetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesSetNumber0

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the time parameters for a time control loop identified by user numbers.
  SUBROUTINE CMISSControlLoopTimesSetNumber1(ProblemUserNumber,ControlLoopIdentifiers,StartTime,StopTime,TimeIncrement,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the time parameters for for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifier.
    REAL(DP), INTENT(IN) :: StartTime !<The start time for the time control loop to set.
    REAL(DP), INTENT(IN) :: StopTime !<The stop time for the time control loop to set.
    REAL(DP), INTENT(IN) :: TimeIncrement !<The time increment for the time control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("CMISSControlLoopTimesSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopTimesSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesSetNumber1

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the time parameters for a time control loop identified by an object.
  SUBROUTINE CMISSControlLoopTimesSetObj(ControlLoop,StartTime,StopTime,TimeIncrement,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(INOUT) :: ControlLoop !<The control loop to set the times for.
    REAL(DP), INTENT(IN) :: StartTime !<The start time for the time control loop to set.
    REAL(DP), INTENT(IN) :: StopTime !<The stop time for the time control loop to set.
    REAL(DP), INTENT(IN) :: TimeIncrement !<The time increment for the time control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTimesSetObj",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_TIMES_SET(ControlLoop%CONTROL_LOOP,StartTime,StopTime,TimeIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTimesSetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesSetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesSetObj")
    CALL CMISS_HANDLE_ERRORS(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesSetObj

  !  
  !================================================================================================================================
  !  
  
  !>Sets/changes the loop type for a control loop identified by user numbers. \todo is this really public???
  SUBROUTINE CMISSControlLoopTypeSetNumber0(ProblemUserNumber,ControlLoopIdentifier,LoopType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the loop type for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    INTEGER(INTG), INTENT(IN) :: LoopType !<The type of control loop to set. \see OPENCMISS_ProblemControlLoopTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopTypeSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,LoopType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopTypeSetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeSetNumber0",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeSetNumber0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeSetNumber0

  !  
  !================================================================================================================================
  !  
  
  !>Sets/changes the loop type for a control loop identified by user numbers. \todo is this really public???
  SUBROUTINE CMISSControlLoopTypeSetNumber1(ProblemUserNumber,ControlLoopIdentifiers,LoopType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to set the loop type for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: LoopType !<The type of control loop to set. \see OPENCMISS_ProblemControlLoopTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSControlLoopTypeSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    IF(ASSOCIATED(PROBLEM)) THEN
      CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
      CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,LoopType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A problem with an user number of "//TRIM(NUMBER_TO_VSTRING(ProblemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSControlLoopTypeSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeSetNumber1

  !  
  !================================================================================================================================
  !  
  
  !>Sets/changes the loop type for a control loop identified by an object. \todo is this really public???
  SUBROUTINE CMISSControlLoopTypeSetObj(ControlLoop,LoopType,Err)
  
    !Argument variables
    TYPE(CMISSControlLoopType), INTENT(INOUT) :: ControlLoop !<The control loop to set the loop type for.
    INTEGER(INTG), INTENT(IN) :: LoopType !<The type of control loop to set. \see OPENCMISS_ProblemControlLoopTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTypeSetObj",Err,ERROR,*999)

    CALL CONTROL_LOOP_TYPE_SET(ControlLoop%CONTROL_LOOP,LoopType,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeSetObj",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeSetObj

!!==================================================================================================================================
!!
!! COORDINATE_ROUTINES
!!
!!==================================================================================================================================
  
  !>Finishes the creation of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemCreateFinishNumber(CoordinateSystemUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to finish creating.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemCreateFinishNumber

  !  
  !================================================================================================================================
  !
  
  !>Finishes the creation of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemCreateFinishObj(CoordinateSystem,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT) :: CoordinateSystem !<The coordinate system to finish creating.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSCoordinateSystemCreateFinishObj",Err,ERROR,*999)
 
    CALL COORDINATE_SYSTEM_CREATE_FINISH(CoordinateSystem%COORDINATE_SYSTEM,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemCreateFinishObj

  !  
  !================================================================================================================================
  !
  
  !>Starts the creation of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemCreateStartNumber(CoordinateSystemUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to start creating.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("CMISSCoordinateSystemCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_CREATE_START(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemCreateStartNumber

  !  
  !================================================================================================================================
  !
  
  !>Starts the creation of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemCreateStartObj(CoordinateSystemUserNumber,CoordinateSystem,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to start creating.
    TYPE(CMISSCoordinateSystemType), INTENT(OUT) :: CoordinateSystem !<On return, the coordinate system that has been created.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSCoordinateSystemCreateStartObj",Err,ERROR,*999)
 
    CALL COORDINATE_SYSTEM_CREATE_START(CoordinateSystemUserNumber,CoordinateSystem%COORDINATE_SYSTEM,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemCreateStartObj

  !  
  !================================================================================================================================
  !

  !>Destroys a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemDestroyNumber(CoordinateSystemUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_DESTROY(COORDINATE_SYSTEM,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDestroyNumber

  !  
  !================================================================================================================================
  !
  
  !>Destroys a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemDestroyObj(CoordinateSystem,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT) :: CoordinateSystem !<The coordinate system to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSCoordinateSysteDestroyObj",Err,ERROR,*999)
 
    CALL COORDINATE_SYSTEM_DESTROY(CoordinateSystem%COORDINATE_SYSTEM,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemDestroyObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDestroyObj",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Returns the dimension of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemDimensionGetNumber(CoordinateSystemUserNumber,CoordinateSystemDimension,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to get the dimension for.
    INTEGER(INTG), INTENT(OUT) :: CoordinateSystemDimension !<On return, the dimension of the coordinate system
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemDimensionGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,CoordinateSystemDimension,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemDimensionGetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDimensionGetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDimensionGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDimensionGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the dimension of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemDimensionGetObj(CoordinateSystem,CoordinateSystemDimension,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(IN) :: CoordinateSystem !<The coordinate system to get the dimension for.
    INTEGER(INTG), INTENT(OUT) :: CoordinateSystemDimension !<On return, the dimension of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemDimensionGetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_DIMENSION_GET(CoordinateSystem,CoordinateSystemDimension,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemDimensionGetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDimensionGetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemDimensionGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDimensionGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the dimension of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemDimensionSetNumber(CoordinateSystemUserNumber,CoordinateSystemDimension,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to set the dimension for.
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemDimension !<The dimension of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemDimensionSetNumber",Err,ERROR,*999)
    
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,CoordinateSystemDimension,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemDimensionSetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDimensionSetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDimensionSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDimensionSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the dimension of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemDimensionSetObj(CoordinateSystem,CoordinateSystemDimension,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT)  :: CoordinateSystem !<The coordinate system to set the dimension for.
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemDimension !<The dimension of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemDimensionSetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_DIMENSION_SET(CoordinateSystem%COORDINATE_SYSTEM,CoordinateSystemDimension,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemDimensionSetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDimensionSetObj",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDimensionSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDimensionSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Returns the focus of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemFocusGetNumber(CoordinateSystemUserNumber,Focus,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to get the focus for.
    REAL(DP), INTENT(OUT) :: Focus !<On return, the focus of the coordinate system
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemFocusGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_FOCUS_GET(COORDINATE_SYSTEM,Focus,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemFocusGetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemFocusGetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemFocusGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemFocusGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the focus of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemFocusGetObj(CoordinateSystem,Focus,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(IN) :: CoordinateSystem !<The coordinate system to get the focus for.
    REAL(DP), INTENT(OUT) :: Focus !<On return, the focus of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemFocusGetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_FOCUS_GET(CoordinateSystem%COORDINATE_SYSTEM,Focus,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemFocusGetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemFocusGetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemFocusGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemFocusGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the focus of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemFocusSetNumber(CoordinateSystemUserNumber,Focus,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to set the focus for.
    REAL(DP), INTENT(IN) :: Focus !<The focus of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemFocusSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_FOCUS_SET(COORDINATE_SYSTEM,Focus,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemFocusSetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemFocusSetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemFocusSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemFocusSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the focus of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemFocusSetObj(CoordinateSystem,Focus,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT) :: CoordinateSystem !<The coordinate system to set the focus for.
    REAL(DP), INTENT(IN) :: Focus !<The focus of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemFocusSetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_FOCUS_SET(CoordinateSystem%COORDINATE_SYSTEM,Focus,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemFocusSetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemFocusSetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemFocusSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemFocusSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Returns the radial interpolation type of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemRadialInterpolationGetNumber(CoordinateSystemUserNumber,RadialInterpolationType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to get the radial interpolation for.
    INTEGER(INTG), INTENT(OUT) :: RadialInterpolationType !<On return, the radial interpolation type of the coordinate system \see OPENCMISS_CoordinateRadialInterpolations
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_GET(COORDINATE_SYSTEM,RadialInterpolationType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemRadialInterpolationGetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemRadialInterpolationGetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemRadialInterpolationGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemRadialInterpolationGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the radial interpolation type of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemRadialInterpolationGetObj(CoordinateSystem,RadialInterpolationType,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT) :: CoordinateSystem !<The coordinate system to get the radial interpolation type for.
    INTEGER(INTG), INTENT(OUT) :: RadialInterpolationType !<On return, the radial interpolation type of the coordinate system. \see OPENCMISS_CoordinateRadialInterpolations
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationGetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_GET(CoordinateSystem%COORDINATE_SYSTEM,RadialInterpolationType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemRadialInterpolationGetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemRadialInterpolationGetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemRadialInterpolationGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemRadialInterpolationGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the radial interpolation type of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemRadialInterpolationSetNumber(CoordinateSystemUserNumber,RadialInterpolationType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to set the radial interpolation for.
    INTEGER(INTG), INTENT(IN) :: RadialInterpolationType !<The radial interpolation type of the coordinate system to set.\see OPENCMISS_CoordinateRadialInterpolations
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_SET(COORDINATE_SYSTEM,RadialInterpolationType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemRadialInterpolationSetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemRadialInterpolationSetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemRadialInterpolationSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemRadialInterpolationSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the radial interpolation type of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemRadialInterpolationSetObj(CoordinateSystem,RadialInterpolationType,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT)  :: CoordinateSystem !<The coordinate system to set the radial interpolation type for.
    INTEGER(INTG), INTENT(IN) :: RadialInterpolationType !<The radial interpolation type of the coordinate system to set. \see OPENCMISS_CoordinateRadialInterpolations
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationSetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_SET(CoordinateSystem%COORDINATE_SYSTEM,RadialInterpolationType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemRadialInterpolationSetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemRadialInterpolationSetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemRadialInterpolationSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemRadialInterpolationSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Returns the type of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemTypeGetNumber(CoordinateSystemUserNumber,CoordinateSystemType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to get the type for.
    INTEGER(INTG), INTENT(OUT) :: CoordinateSystemType !<On return, the type of the coordinate system. \see OPENCMISS_CoordinateSystemTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_TYPE_GET(COORDINATE_SYSTEM,CoordinateSystemType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the type of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemTypeGetObj(CoordinateSystem,CoordinateSystemType,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(IN) :: CoordinateSystem !<The coordinate system to get the type for.
    INTEGER(INTG), INTENT(OUT) :: CoordinateSystemType !<On return, the type of the coordinate system. \see OPENCMISS_CoordinateSystemTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemTypeGetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_TYPE_GET(CoordinateSystem%COORDINATE_SYSTEM,CoordinateSystemType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeGetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the type of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemTypeSetNumber(CoordinateSystemUserNumber,CoordinateSystemType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to set the type for.
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemType !<The type of the coordinate system to set. \see OPENCMISS_CoordinateSystemTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemTypeSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_TYPE_SET(COORDINATE_SYSTEM,CoordinateSystemType,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemTypeSetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeSetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemTypeSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the type of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemTypeSetObj(CoordinateSystem,CoordinateSystemType,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT) :: CoordinateSystem !<The coordinate system to set the type for.
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemType !<The type of the coordinate system to set. \see OPENCMISS_CoordinateSystemTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemTypeSetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_TYPE_SET(CoordinateSystem%COORDINATE_SYSTEM,CoordinateSystemType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeSetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeSetObj
  
  !  
  !================================================================================================================================
  !  
 
  !>Returns the origin of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemOriginGetNumber(CoordinateSystemUserNumber,Origin,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to get the origin for.
    REAL(DP), INTENT(OUT) :: Origin(:) !<On return, the orign of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemOriginGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_ORIGIN_GET(COORDINATE_SYSTEM,Origin,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemOriginGetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOriginGetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemOriginGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOriginGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the origin of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemOriginGetObj(CoordinateSystem,Origin,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(IN) :: CoordinateSystem !<The coordinate system to get the origin for.
    REAL(DP), INTENT(OUT) :: Origin(:) !<On return, the origin of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOriginGetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIGIN_GET(CoordinateSystem%COORDINATE_SYSTEM,Origin,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOriginGetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOriginGetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOriginGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOriginGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the origin of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemOriginSetNumber(CoordinateSystemUserNumber,Origin,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to set the origin for.
    REAL(DP), INTENT(IN) :: Origin(:) !<The orign of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemOriginSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_ORIGIN_SET(COORDINATE_SYSTEM,Origin,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemOriginSetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOriginSetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemOriginSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOriginSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the origin of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemOriginSetObj(CoordinateSystem,Origin,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(IN) :: CoordinateSystem !<The coordinate system to set the origin for.
    REAL(DP), INTENT(IN) :: Origin(:) !<The origin of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOriginSetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIGIN_SET(CoordinateSystem%COORDINATE_SYSTEM,Origin,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOriginSetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOriginSetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOriginSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOriginSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Returns the orientation of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemOrientationGetNumber(CoordinateSystemUserNumber,Orientation,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to get the orientation for.
    REAL(DP), INTENT(OUT) :: Orientation(:,:) !<On return, the orientation of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemOrientationGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_ORIENTATION_GET(COORDINATE_SYSTEM,ORIENTATION,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemOrientationGetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOrientationGetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemOrientationGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOrientationGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the orientation of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemOrientationGetObj(CoordinateSystem,Orientation,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(IN) :: CoordinateSystem !<The coordinate system to get the orientation for.
    REAL(DP), INTENT(OUT) :: Orientation(:,:) !<On return, the orientation of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOrientationGetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIENTATION_GET(CoordinateSystem%COORDINATE_SYSTEM,Orientation,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOrientationGetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOrientationGetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOrientationGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOrientationGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the orientation of a coordinate system identified by a user number.
  SUBROUTINE CMISSCoordinateSystemOrientationSetNumber(CoordinateSystemUserNumber,Orientation,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to set the orientation for.
    REAL(DP), INTENT(IN) :: Orientation(:,:) !<The orientation of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSCoordinateSystemOrientationSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL COORDINATE_SYSTEM_ORIENTATION_SET(COORDINATE_SYSTEM,ORIENTATION,Err,ERROR,*999)
    ELSE
      LOCAL_ERROR="A coordinate system with an user number of "// &
        & TRIM(NUMBER_TO_VSTRING(CoordinateSystemUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSCoordinateSystemOrientationSetNumber")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOrientationSetNumber",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemOrientationSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOrientationSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the orientation of a coordinate system identified by an object.
  SUBROUTINE CMISSCoordinateSystemOrientationSetObj(CoordinateSystem,Orientation,Err)
  
    !Argument variables
    TYPE(CMISSCoordinateSystemType), INTENT(INOUT) :: CoordinateSystem !<The coordinate system to set the orientation for.
    REAL(DP), INTENT(IN) :: Orientation(:,:) !<The orientation of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOrientationSetObj",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIENTATION_SET(CoordinateSystem%COORDINATE_SYSTEM,Orientation,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOrientationSetObj")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOrientationSetObj",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOrientationSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOrientationSetObj

!!==================================================================================================================================
!!
!! EQUATIONS_ROUTINES
!!
!!==================================================================================================================================

  !>Destroys equations for equations identified by a user number.
  SUBROUTINE CMISSEquationsDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to destroy.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_DESTROY(EQUATIONS,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsDestroyNumber

  !  
  !================================================================================================================================
  !  
 
  !>Destroy equations for equations identified by an object.
  SUBROUTINE CMISSEquationsDestroyObj(Equations,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(INOUT) :: Equations !<The equations to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_DESTROY(Equations%EQUATIONS,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsDestroyObj
  
  !  
  !================================================================================================================================
  !  
 
  !>Gets the linearity type for equations identified by a user number.
  SUBROUTINE CMISSEquationsLinearityTypeGetNumber(RegionUserNumber,EquationsSetUserNumber,LinearityType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to get the linearity type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the linearity type for.
    INTEGER(INTG), INTENT(OUT) :: LinearityType !<On return, the linearity type of the equations \see OPENCMISS_EquationsLinearityTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsLinearityTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_LINEARITY_TYPE_GET(EQUATIONS,LinearityErr,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsLinearityTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsLinearityTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsLinearityTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsLinearityTypeGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Gets the linearity type for equations identified by an object.
  SUBROUTINE CMISSEquationsLinearityTypeGetObj(Equations,LinearityType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(IN) :: Equations !<The equations to get the linearity type for.
    INTEGER(INTG), INTENT(OUT) :: LinearityType !<On return, the linearity type of the equations \see OPENCMISS_EquationsLinearityTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsLinearityTypeGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_LINEARITY_TYPE_GET(Equations%EQUATIONS,LinearityType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsLinearityTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsLinearityTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsLinearityTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsLinearityTypeGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Gets the lumping type for equations identified by a user number.
  SUBROUTINE CMISSEquationsLumpingTypeGetNumber(RegionUserNumber,EquationsSetUserNumber,LumpingType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to get the lumping type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the lumping type for.
    INTEGER(INTG), INTENT(OUT) :: LumpingType !<On return, the lumping type of the equations \see OPENCMISS_EquationsLumpingTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsLumpingTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_LUMPING_TYPE_GET(EQUATIONS,LumpingType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsLumpingTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsLumpingTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsLumpingTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsLumpingTypeGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Gets the lumping type for equations identified by an object.
  SUBROUTINE CMISSEquationsLumpingTypeGetObj(Equations,LumpingType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(IN) :: Equations !<The equations to get the lumping type for.
    INTEGER(INTG), INTENT(OUT) :: LumpingType !<On return, the lumping type of the equations \see OPENCMISS_EquationsLumpingTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsLumpingTypeGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_LUMPING_TYPE_GET(Equations%EQUATIONS,LumpingType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsLumpingTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsLumpingTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsLumpingTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsLumpingTypeGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the lumping type for equations identified by a user number.
  SUBROUTINE CMISSEquationsLumpingTypeSetNumber(RegionUserNumber,EquationsSetUserNumber,LumpingType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to set the lumping type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the lumping type for.
    INTEGER(INTG), INTENT(IN) :: LumpingType !<The lumping type of the equations to set\see OPENCMISS_EquationsLumpingTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsLumpingTypeSetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_LUMPING_TYPE_SET(EQUATIONS,LumpingType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsLumpingTypeSetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsLumpingTypeSetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsLumpingTypeSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsLumpingTypeSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the lumping type for equations identified by an object.
  SUBROUTINE CMISSEquationsLumpingTypeSetObj(Equations,LumpingType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(INOUT) :: Equations !<The equations to set the lumping type for.
    INTEGER(INTG), INTENT(IN) :: LumpingType !<The lumping type of the equations to set\see OPENCMISS_EquationsLumpingTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsLumpingTypeSetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_LUMPING_TYPE_SET(Equations%EQUATIONS,LumpingType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsLumpingTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsLumpingTypeSetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsLumpingTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsLumpingTypeSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Gets the output type for equations identified by a user number.
  SUBROUTINE CMISSEquationsOutputTypeGetNumber(RegionUserNumber,EquationsSetUserNumber,OutputType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to get the output type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the output type for.
    INTEGER(INTG), INTENT(OUT) :: OutputType !<On return, the output type of the equations \see OPENCMISS_EquationsOutputTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsOutputTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_OUTPUT_TYPE_GET(EQUATIONS,OutputType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsOutputTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsOutputTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsOutputTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsOutputTypeGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Gets the output type for equations identified by an object.
  SUBROUTINE CMISSEquationsOutputTypeGetObj(Equations,OutputType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(IN) :: Equations !<The equations to get the output type for.
    INTEGER(INTG), INTENT(OUT) :: OutputType !<On return, the output type of the equations \see OPENCMISS_EquationsOutputTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsOutputTypeGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_OUTPUT_TYPE_GET(Equations%EQUATIONS,OutputType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsOutputTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsOutputTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsOutputTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsOutputTypeGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the output type for equations identified by a user number.
  SUBROUTINE CMISSEquationsOutputTypeSetNumber(RegionUserNumber,EquationsSetUserNumber,OutputType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to set the output type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the output type for.
    INTEGER(INTG), INTENT(IN) :: OutputType !<The output type of the equations to set \see OPENCMISS_EquationsOutputTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsOutputTypeSetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,OutputType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsOutputTypeSetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsOutputTypeSetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsOutputTypeSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsOutputTypeSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the output type for equations identified by an object.
  SUBROUTINE CMISSEquationsOutputTypeSetObj(Equations,OutputType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(INOUT) :: Equations !<The equations to set the output type for.
    INTEGER(INTG), INTENT(IN) :: OutputType !<The output type of the equations to set \see OPENCMISS_EquationsOutputTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsOutputTypeSetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_OUTPUT_TYPE_SET(Equations%EQUATIONS,OutputType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsOutputTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsOutputTypeSetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsOutputTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsOutputTypeSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Gets the sparsity type for equations identified by a user number.
  SUBROUTINE CMISSEquationsSparsityTypeGetNumber(RegionUserNumber,EquationsSetUserNumber,SparsityType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to get the sparsity type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the sparsity type for.
    INTEGER(INTG), INTENT(OUT) :: SparsityType !<On return, the sparsity type of the equations \see OPENCMISS_EquationsSparsityTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSparsityTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_SPARSITY_TYPE_GET(EQUATIONS,SparsityType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSparsityTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSparsityTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSparsityTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSparsityTypeGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Gets the sparsity type for equations identified by an object.
  SUBROUTINE CMISSEquationsSparsityTypeGetObj(Equations,SparsityType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(IN) :: Equations !<The equations to get the sparsity type for.
    INTEGER(INTG), INTENT(OUT) :: SparsityType !<On return, the sparsity type of the equations \see OPENCMISS_EquationsSparsityTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSparsityTypeGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SPARSITY_TYPE_GET(Equations%EQUATIONS,SparsityType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSparsityTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSparsityTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSparsityTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSparsityTypeGetObj

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the sparsity type for equations identified by a user number.
  SUBROUTINE CMISSEquationsSparsityTypeSetNumber(RegionUserNumber,EquationsSetUserNumber,SparsityType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to set the sparsity type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the sparsity type for.
    INTEGER(INTG), INTENT(IN) :: SparsityType !<The sparsity type of the equations to set \see OPENCMISS_EquationsSparsityTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSparsityTypeSetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,SparsityType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSparsityTypeSetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSparsityTypeSetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSparsityTypeSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSparsityTypeSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the sparsity type for equations identified by an object.
  SUBROUTINE CMISSEquationsSparsityTypeSetObj(Equations,SparsityType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(INOUT) :: Equations !<The equations to set the sparsity type for.
    INTEGER(INTG), INTENT(IN) :: SparsityType !<The sparsity type of the equations to set \see OPENCMISS_EquationsSparsityTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSparsityTypeSetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SPARSITY_TYPE_SET(Equations%EQUATIONS,SparsityType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSparsityTypeSetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSparsityTypeSetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSparsityTypeSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSparsityTypeSetObj

  !  
  !================================================================================================================================
  !  
 
  !>Gets the time dependence type for equations identified by a user number.
  SUBROUTINE CMISSEquationsTimeDependenceTypeGetNumber(RegionUserNumber,EquationsSetUserNumber,TimeDependenceType,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations to get the time dependence type for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the time dependence type for.
    INTEGER(INTG), INTENT(OUT) :: TimeDependenceType !<On return, the time dependence type of the equations \see OPENCMISS_EquationsTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsTimeDependenceTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
        CALL EQUATIONS_TIME_DEPENDENCE_TYPE_GET(EQUATIONS,TimeDependenceType,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsTimeDependenceTypeGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsTimeDependenceTypeGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsTimeDependenceTypeGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsTimeDependenceTypeGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Gets the time dependence type for equations identified by an object.
  SUBROUTINE CMISSEquationsTimeDependenceTypeGetObj(Equations,TimeDependenceType,Err)
  
    !Argument variables
    TYPE(CMISSEquationsType), INTENT(IN) :: Equations !<The equations to get the time dependence type for.
    INTEGER(INTG), INTENT(OUT) :: TimeDependenceType !<On return, the time dependence type of the equations \see OPENCMISS_EquationsTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsTimeDependenceTypeGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_TIME_DEPENDENCE_TYPE_GET(Equations%EQUATIONS,TimeDependenceType,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsTimeDependenceTypeGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsTimeDependenceTypeGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsTimeDependenceTypeGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsTimeDependenceTypeGetObj


!!==================================================================================================================================
!!
!! EQUATIONS_SET_ROUTINES
!!
!!==================================================================================================================================

  !>Finish the creation of a analytic solution for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetAnalyticCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetAnalyticCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_ANALYTIC_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetAnalyticCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetAnalyticCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetAnalyticCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetAnalyticCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of a analytic solution for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetAnalyticCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetAnalyticCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_ANALYTIC_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetAnalyticCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetAnalyticCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetAnalyticCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetAnalyticCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of a analytic solution for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetAnalyticCreateStartNumber(RegionUserNumber,EquationsSetUserNumber,AnalyticFunctionType, &
    & AnalyticFieldUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of.
    INTEGER(INTG), INTENT(IN) :: AnalyticFunctionType !<The analytic function type to use. \see OPENCMISS_EquationsSetAnalyticFunctionTypes
    INTEGER(INTG), INTENT(IN) :: AnalyticFieldUserNumber !<The user number of the field for the analytic function
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetAnalyticCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(ANALYTIC_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL FIELD_USER_NUMBER_FIND(AnalyticFieldUserNumber,REGION,ANALYTIC_FIELD,ERR,ERROR,*999)
        CALL EQUATIONS_SET_ANALYTIC_CREATE_START(EQUATIONS_SET,AnalyticFunctionType,AnalyticFieldUserNumber,ANALYTIC_FIELD, &
          & Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetAnalyticCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetAnalyticCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetAnalyticCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetAnalyticCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of an analytic solution for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetAnalyticCreateStartObj(EquationsSet,AnalyticFunctionType,AnalyticFieldUserNumber,AnalyticField,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<The equations set to start the analytic creation on.
    INTEGER(INTG), INTENT(IN) :: AnalyticFunctionType !<The analytic function type to use. \see OPENCMISS_EquationsSetAnalyticFunctionTypes
    INTEGER(INTG), INTENT(IN) :: AnalyticFieldUserNumber !<The user number of the field for the analytic function
    TYPE(CMISSFieldType), INTENT(INOUT) :: AnalyticField !<If associated on entry, the user created analytic field which has the same user number as the specified analytic field user number. If not associated on entry, on return, the created analytic field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetAnalyticCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_ANALYTIC_CREATE_START(EquationsSet%EQUATIONS_SET,AnalyticFunctionType,AnalyticFieldUserNumber, &
      & AnalyticField%FIELD,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetAnalyticCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetAnalyticCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetAnalyticCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetAnalyticCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy the analytic solution for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetAnalyticDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destroy.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
     TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetAnalyticDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(ANALYTIC_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_ANALYTIC_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetAnalyticDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetAnalyticDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetAnalyticDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetAnalyticDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy the analytic solution for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetAnalyticDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy the analytic for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetAnalyticDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_ANALYTIC_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetAnalyticDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetAnalyticDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetAnalyticDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetAnalyticDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Set boundary conditions for an equation set according to the analytic equations for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsAnalyticNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to set the analytic boundary conditions.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the analytic boundary conditions.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetBoundaryConditionsAnalyticNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(ANALYTIC_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetBoundaryConditionsAnalyticNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsAnalyticNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsAnalyticNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsAnalyticNumber

  !  
  !================================================================================================================================
  !  

  !>Set boundary conditions for an equation set according to the analytic equations for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsAnalyticObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to set the analytic boundary conditions.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsAnalyticObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetBoundaryConditionsAnalyticObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsAnalyticObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsAnalyticObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsAnalyticObj

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of boundary conditions for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the boundary conditions to finish.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of a boundary conditions for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish the creation of boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of boundary conditions for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateStartNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the boundary conditions to start the creation of.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to start the creation of boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetBoundaryConditionsCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of boundary conditions for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateStartObj(EquationsSet,BoundaryConditions,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<The equations set to start the creation of boundary conditions on.
    TYPE(CMISSBoundaryConditionsType), INTENT(INOUT) :: BoundaryConditions !<On return, the created boundary conditions.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EquationsSet%EQUATIONS_SET,BoundaryConditions%BOUNDARY_CONDITIONS, &
      & Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy the boundary conditions for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destory the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetBoundaryConditionsDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetBoundaryConditionsDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy the boundary conditions for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetBoundaryConditionsDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy the boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetBoundaryConditionsDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsDestroyObj

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetCreateStartNumber(EquationsSetUserNumber,RegionUserNumber,GeomFibreFieldUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to be created.
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region to start the creation of an equations set on.
    INTEGER(INTG), INTENT(IN) :: GeomFibreFieldUserNumber !<The user number of the Geometric/Fibre field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: GEOM_FIBRE_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(GEOM_FIBRE_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL FIELD_USER_NUMBER_FIND(GeomFibreFieldUserNumber,REGION,GEOM_FIBRE_FIELD,Err,ERROR,*999)
      IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
        CALL EQUATIONS_SET_CREATE_START(EquationsSetUserNumber,REGION,GEOM_FIBRE_FIELD,EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="A field with an user number of "//TRIM(NUMBER_TO_VSTRING(GeomFibreFieldUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetCreateStartObj(EquationsSetUserNumber,Region,GeomFibreField,EquationsSet,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to be created.
    TYPE(CMISSRegionType), INTENT(IN) :: Region !<The region to create the equations set on.
    TYPE(CMISSFieldType), INTENT(IN) :: GeomFibreField !<The Geometric/Fibre field for the creation of the equations set.
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<On return, the created equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_CREATE_START(EquationsSetUserNumber,Region%REGION,GeomFibreField%FIELD,EquationsSet%EQUATIONS_SET, &
      & Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destory.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
     TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Finish the creation of dependent variables for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetDependentCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish the creation of dependent variables for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of dependent variables for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetDependentCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetDependentCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDependentCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDependentCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDependentCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of dependent variables for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetDependentCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish the creation of dependent variables for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetDependentCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetDependentCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDependentCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDependentCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDependentCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of dependent variables for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetDependentCreateStartNumber(RegionUserNumber,EquationsSetUserNumber,DependentFieldUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to start the creation of dependent variables for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to start the creation of dependent variables for.
    INTEGER(INTG), INTENT(IN) :: DependentFieldUserNumber !<The user number of the dependent field.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetDependentCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(DEPENDENT_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL FIELD_USER_NUMBER_FIND(DependentFieldUserNumber,REGION,DEPENDENT_FIELD,ERR,ERROR,*999)
        CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,DependentFieldUserNumber,DEPENDENT_FIELD,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetDependentCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDependentCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDependentCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDependentCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of dependent variables for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetDependentCreateStartObj(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<The equations set to start the creation of dependent variables on.
    INTEGER(INTG), INTENT(IN) :: DependentFieldUserNumber !<The user number of the dependent field.
    TYPE(CMISSFieldType), INTENT(INOUT) :: DependentField !<If associated on entry, the user created dependent field which has the same user number as the specified dependent field user number. If not associated on entry, on return, the created dependent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetDependentCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EquationsSet%EQUATIONS_SET,DependentFieldUserNumber,DependentField%FIELD, &
      & Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetDependentCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDependentCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDependentCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDependentCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy the dependent variables for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetDependentDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destroy the dependent variables for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the dependent variables for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetDependentDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_DEPENDENT_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetDependentDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDependentDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDependentDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDependentDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy the dependent variables for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetDependentDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy the dependent variables for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetDependentDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_DEPENDENT_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetDependentDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetDependentDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetDependentDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetDependentDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Finish the creation of equations for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetEquationsCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish the creation of equations for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of equations for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetEquationsCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetEquationsCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetEquationsCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetEquationsCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetEquationsCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of equations for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetEquationsCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish the creation of equations for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetEquationsCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetEquationsCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetEquationsCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetEquationsCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetEquationsCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of equations for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetEquationsCreateStartNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to start the creation of equations for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to start the creation of equations for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetEquationsCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL FIELD_USER_NUMBER_FIND(DependentFieldUserNumber,REGION,DEPENDENT_FIELD,ERR,ERROR,*999)
        CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetEquationsCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetEquationsCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetEquationsCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetEquationsCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of equations for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetEquationsCreateStartObj(EquationsSet,Equations,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<The equations set to start the creation of equations on.
    TYPE(CMISSEquationsType), INTENT(INOUT) :: Equations !<On return, the created equations.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetEquationsCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EquationsSet%EQUATIONS_SET,Equations%EQUATIONS,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetEquationsCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetEquationsCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetEquationsCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetEquationsCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy the equations for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetEquationsDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destroy the equations for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetEquationsDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_EQUATIONS_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetEquationsDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetEquationsDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetEquationsDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetEquationsDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy the equations for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetEquationsDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetEquationsDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_EQUATIONS_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetEquationsDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetEquationsDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetEquationsDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetEquationsDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Finish the creation of materials for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetMaterialsCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish the creation of materials for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of materials for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetMaterialsCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetMaterialsCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetMaterialsCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetMaterialsCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetMaterialsCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of materials for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetMaterialsCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish the creation of materials for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetMaterialsCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_MATERIALS_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetMaterialsCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetMaterialsCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetMaterialsCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetMaterialsCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of materials for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetMaterialsCreateStartNumber(RegionUserNumber,EquationsSetUserNumber,MaterialsFieldUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to start the creation of materials for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to start the creation of materials for.
    INTEGER(INTG), INTENT(IN) :: MaterialsFieldUserNumber !<The user number of the materials field.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetMaterialsCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(MATERIALS_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL FIELD_USER_NUMBER_FIND(MaterialsFieldUserNumber,REGION,MATERIALS_FIELD,ERR,ERROR,*999)
        CALL EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,MaterialsFieldUserNumber,MATERIALS_FIELD,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetMaterialsCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetMaterialsCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetMaterialsCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetMaterialsCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of materials for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetMaterialsCreateStartObj(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<The equations set to start the creation of materials on.
    INTEGER(INTG), INTENT(IN) :: MaterialsFieldUserNumber !<The user number of the materials field.
    TYPE(CMISSFieldType), INTENT(INOUT) :: MaterialsField !<If associated on entry, the user created materials field which has the same user number as the specified materials field user number. If not associated on entry, on return, the created materials field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetMaterialsCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_MATERIALS_CREATE_START(EquationsSet%EQUATIONS_SET,MaterialsFieldUserNumber,MaterialsField%FIELD, &
      & Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetMaterialsCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetMaterialsCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetMaterialsCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetMaterialsCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy the materials for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetMaterialsDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destroy the materials for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the materials for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetMaterialsDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_MATERIALS_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetMaterialsDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetMaterialsDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetMaterialsDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetMaterialsDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy the materials for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetMaterialsDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy the materials for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetMaterialsDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_MATERIALS_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetMaterialsDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetMaterialsDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetMaterialsDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetMaterialsDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Returns the solution method for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSolutionMethodGetNumber(RegionUserNumber,EquationsSetUserNumber,SolutionMethod,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to get the solution method for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the solution method for.
    INTEGER(INTG), INTENT(OUT) :: SolutionMethod !<On return, the solution method. \see OPENCMISS_EquationsSetSolutionMethods
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetSolutionMethodGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SOLUTION_METHOD_GET(EQUATIONS_SET,SolutionMethod,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSolutionMethodGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSolutionMethodGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSolutionMethodGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSolutionMethodGetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the solution method for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSolutionMethodGetObj(EquationsSet,SolutionMethod,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to get the solution method for.
    INTEGER(INTG), INTENT(OUT) :: SolutionMethod !<On Return, the solution method. \see OPENCMISS_EquationsSetSolutionMethods
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetSolutionMethodGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SOLUTION_METHOD_GET(EquationsSet%EQUATIONS_SET,SolutionMethod,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetSolutionMethodGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSolutionMethodGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSolutionMethodGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSolutionMethodGetObj

  !  
  !================================================================================================================================
  !  

  !>Sets/changes the solution method for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSolutionMethodSetNumber(RegionUserNumber,EquationsSetUserNumber,SolutionMethod,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to set the solution method for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the solution method for.
    INTEGER(INTG), INTENT(IN) :: SolutionMethod !<The solution method to set. \see OPENCMISS_EquationsSetSolutionMethods
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetSolutionMethodSetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SolutionMethod,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSolutionMethodSetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSolutionMethodSetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSolutionMethodSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSolutionMethodSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the solution method for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSolutionMethodSetObj(EquationsSet,SolutionMethod,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to set the solution method for.
    INTEGER(INTG), INTENT(IN) :: SolutionMethod !<The solution method to set. \see OPENCMISS_EquationsSetSolutionMethods
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetSolutionMethodSetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SOLUTION_METHOD_SET(EquationsSet%EQUATIONS_SET,SolutionMethod,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetSolutionMethodSetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSolutionMethodSetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSolutionMethodSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSolutionMethodSetObj

  !  
  !================================================================================================================================
  !  

  !>Finish the creation of a source for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSourceCreateFinishNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to finish the creation of a source for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to finish the creation of a source for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetSourceCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SOURCE_CREATE_FINISH(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSourceCreateFinishNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSourceCreateFinishNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSourceCreateFinishNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSourceCreateFinishNumber

  !  
  !================================================================================================================================
  !  
 
  !>Finish the creation of a source for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSourceCreateFinishObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to finish the creation of a source for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetSourceCreateFinishObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SOURCE_CREATE_FINISH(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetSourceCreateFinishObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSourceCreateFinishObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSourceCreateFinishObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSourceCreateFinishObj

  !  
  !================================================================================================================================
  !  

  !>Start the creation of a source for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSourceCreateStartNumber(RegionUserNumber,EquationsSetUserNumber,SourceFieldUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to start the creation of a source for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to start the creation of a source for.
    INTEGER(INTG), INTENT(IN) :: SourceFieldUserNumber !<The user number of the source field.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetSourceCreateStartNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(SOURCE_FIELD)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL FIELD_USER_NUMBER_FIND(SourceFieldUserNumber,REGION,SOURCE_FIELD,ERR,ERROR,*999)
        CALL EQUATIONS_SET_SOURCE_CREATE_START(EQUATIONS_SET,SourceFieldUserNumber,SOURCE_FIELD,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSourceCreateStartNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSourceCreateStartNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSourceCreateStartNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSourceCreateStartNumber

  !  
  !================================================================================================================================
  !  

  !>Start the creation of a source for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSourceCreateStartObj(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(OUT) :: EquationsSet !<The equations set to start the creation of a source on.
    INTEGER(INTG), INTENT(IN) :: SourceFieldUserNumber !<The user number of the source field.
    TYPE(CMISSFieldType), INTENT(INOUT) :: SourceField !<If associated on entry, the user created source field which has the same user number as the specified source field user number. If not associated on entry, on return, the created source field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetSourceCreateStartObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SOURCE_CREATE_START(EquationsSet%EQUATIONS_SET,SourceFieldUserNumber,SourceField%FIELD,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetSourceCreateStartObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSourceCreateStartObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSourceCreateStartObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSourceCreateStartObj

  !  
  !================================================================================================================================
  !  

  !>Destroy the source for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSourceDestroyNumber(RegionUserNumber,EquationsSetUserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to destroy the source for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to destroy the source for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CMISSEquationsSetSourceDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SOURCE_DESTROY(EQUATIONS_SET,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSourceDestroyNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSourceDestroyNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSourceDestroyNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSourceDestroyNumber

  !  
  !================================================================================================================================
  !  

  !>Destroy the source for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSourceDestroyObj(EquationsSet,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to destroy the source for.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSEquationsSetSourceDestroyObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SOURCE_DESTROY(EquationsSet%EQUATIONS_SET,Err,ERROR,*999)
    
    CALL EXITS("CMISSEquationsSetSourceDestroyObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSourceDestroyObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSourceDestroyObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSourceDestroyObj

  !  
  !================================================================================================================================
  !  

  !>Returns the equations set specification i.e., equations set class, type and subtype for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSpecificationGetNumber(RegionUserNumber,EquationsSetUserNumber,EquationsSetClass, &
    & EquationsSetType,EquationsSetSubtype,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to get the specification for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the specification for.
    INTEGER(INTG), INTENT(OUT) :: EquationsSetClass !<On return, the equations set class. \see OPENCMISS_EquationsSetClasses
    INTEGER(INTG), INTENT(OUT) :: EquationsSetType !<On return, the equations set type. \see OPENCMISS_EquationsSetTypes
    INTEGER(INTG), INTENT(OUT) :: EquationsSetSubtype !<On return, the equations set subtype. \see OPENCMISS_EquationsSetSubtypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetSpecificationGetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SPECIFICATION_GET(EQUATIONS_SET,EquationsSetClass,EquationsSetType,EquationsSetSubtype,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSpecificationGetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSpecificationGetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSpecificationGetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSpecificationSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Returns the equations set specification i.e., equations set class, type and subtype for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSpecificationGetObj(EquationsSet,EquationsSetClass,EquationsSetType,EquationsSetSubtype,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(IN) :: EquationsSet !<The equations set to get the specification for.
    INTEGER(INTG), INTENT(OUT) :: EquationsSetClass !<On return, the equations set class. \see OPENCMISS_EquationsSetClasses
    INTEGER(INTG), INTENT(OUT) :: EquationsSetType !<On return, the equations set type. \see OPENCMISS_EquationsSetTypes
    INTEGER(INTG), INTENT(OUT) :: EquationsSetSubtype !<On return, the equations set subtype. \see OPENCMISS_EquationsSetSubtypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetSpecificationGetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SPECIFICATION_GET(EquationsSet%EQUATIONS_SET,EquationsSetClass,EquationsSetType,EquationsSetSubtype, &
      & Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetSpecificationGetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSpecificationGetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSpecificationGetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSpecificationGetObj

  !  
  !================================================================================================================================
  !  

  !>Sets/changes the equations set specification i.e., equations set class, type and subtype for an equations set identified by a user number.
  SUBROUTINE CMISSEquationsSetSpecificationSetNumber(RegionUserNumber,EquationsSetUserNumber,EquationsSetClass, &
    & EquationsSetType,EquationsSetSubtype,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the Region containing the equations set to set the specification for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to set the specification for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetClass !<The equations set class to set. \see OPENCMISS_EquationsSetClasses
    INTEGER(INTG), INTENT(IN) :: EquationsSetType !<The equations set type to set. \see OPENCMISS_EquationsSetTypes
    INTEGER(INTG), INTENT(IN) :: EquationsSetSubtype !<The equations set subtype to set. \see OPENCMISS_EquationsSetSubtypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CMISSEquationsSetSpecificationSetNumber",Err,ERROR,*999)
 
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(REGION,EquationsSetUserNumber,Err,EQUATIONS_SET,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EquationsSetClass,EquationsSetType,EquationsSetSubtype,Err,ERROR,*999)
      ELSE
        LOCAL_ERROR="An equations set with an user number of "//TRIM(NUMBER_TO_VSTRING(EquationsSetUserNumber,"*",Err,ERROR))// &
          & " does not exist on region number "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="A region with an user number of "//TRIM(NUMBER_TO_VSTRING(RegionUserNumber,"*",Err,ERROR))//" does not exist."
      CALL FLAG_ERROR(LOCAL_ERROR,Err,ERROR,*999)
    ENDIF

    CALL EXITS("CMISSEquationsSetSpecificationSetNumber")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSpecificationSetNumber",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSpecificationSetNumber")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSpecificationSetNumber

  !  
  !================================================================================================================================
  !  
 
  !>Sets/changes the equations set specification i.e., equations set class, type and subtype for an equations set identified by an object.
  SUBROUTINE CMISSEquationsSetSpecificationSetObj(EquationsSet,EquationsSetClass,EquationsSetType,EquationsSetSubtype,Err)
  
    !Argument variables
    TYPE(CMISSEquationsSetType), INTENT(INOUT) :: EquationsSet !<The equations set to set the specification for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetClass !<The equations set class to set. \see OPENCMISS_EquationsSetClasses
    INTEGER(INTG), INTENT(IN) :: EquationsSetType !<The equations set type to set. \see OPENCMISS_EquationsSetTypes
    INTEGER(INTG), INTENT(IN) :: EquationsSetSubtype !<The equations set subtype to set. \see OPENCMISS_EquationsSetSubtypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSEquationsSetSpecificationSetObj",Err,ERROR,*999)
 
    CALL EQUATIONS_SET_SPECIFICATION_SET(EquationsSet%EQUATIONS_SET,EquationsSetClass,EquationsSetType,EquationsSetSubtype, &
      & Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetSpecificationSetObj")
    RETURN
999 CALL ERRORS("CMISSEquationsSetSpecificationSetObj",Err,ERROR)
    CALL EXITS("CMISSEquationsSetSpecificationSetObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetSpecificationSetObj

!!==================================================================================================================================
!!
!! FIELD_ROUTINES
!!
!!==================================================================================================================================


!!==================================================================================================================================
!!
!! FIELD_IO_ROUTINES
!!
!!==================================================================================================================================

  !>Export element information for fields set identified by an object. \todo number method
  SUBROUTINE CMISSFieldIOElementsExportObj(Fields,FileName,Method,Err)
  
    !Argument variables
    TYPE(CMISSFieldsType), INTENT(INOUT) :: Fields !<The fields to export the elements for.
    CHARACTER(LEN=*), INTENT(IN) :: FileName !<The file name to export the elements to
    CHARACTER(LEN=*), INTENT(IN):: Method !<The export method to use.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSFieldIOElementsExportObj",Err,ERROR,*999)
 
    CALL FIELD_IO_NODES_EXPORT(Fields%FIELDS,FileName,Method,Err,ERROR,*999)

    CALL EXITS("CMISSFieldIOElementsExportObj")
    RETURN
999 CALL ERRORS("CMISSFieldIOElementsExportObj",Err,ERROR)
    CALL EXITS("CMISSFieldIOElementsExportObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFieldIOElementsExportObj

  !  
  !================================================================================================================================
  !  

  !>Export nodal information for fields set identified by an object. \todo number method
  SUBROUTINE CMISSFieldIONodesExportObj(Fields,FileName,Method,Err)
  
    !Argument variables
    TYPE(CMISSFieldsType), INTENT(INOUT) :: Fields !<The fields to export the nodes for.
    CHARACTER(LEN=*), INTENT(IN) :: FileName !<The file name to export the nodes to
    CHARACTER(LEN=*), INTENT(IN):: Method !<The export method to use.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSFieldIONodesExportObj",Err,ERROR,*999)
 
    CALL FIELD_IO_NODES_EXPORT(Fields%FIELDS,FileName,Method,Err,ERROR,*999)

    CALL EXITS("CMISSFieldIONodesExportObj")
    RETURN
999 CALL ERRORS("CMISSFieldIONodesExportObj",Err,ERROR)
    CALL EXITS("CMISSFieldIONodesExportObj")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSFieldIONodesExportObj

  !  
  !================================================================================================================================
  !  
 

 
END MODULE OPENCMISS
