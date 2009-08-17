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
  USE BOUNDARY_CONDITION_ROUTINES
  USE CMISS_ROUTINES
  USE COOORDINATE_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
  USE PROBLEM_CONSTANTS
  
  IMPLICIT NONE

  PRIVATE
  
  !Module parameters
  
  !Module types

  !Module variables

  TYPE(VARYING_STRING) :: ERROR

  !Interfaces

  INTERFACE CMISSInitialise
    MODULE PROCEDURE CMISSInitialiseNumber
    MOUDLE PROCEDURE CMISSInitialisePtr
  END INTERFACE !CMISSInitialise
  
  PUBLIC CMISSFinalise,CMISSInitialise

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
    MODULE PROCEDURE CMISSAnalyticAnalysisOutputPtr
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
  INTEGER(INTG), PARAMETER :: CMISSBasisCubicHermiteInterpolation = BASIS_CUBIC_HERMITE_INTERPOLATION=4 !<Cubic Hermite interpolation specification \see OPENCMISS_BasisInterpolationSpecifications,OPENCMISS
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
    MODULE PROCEDURE CMISSBasisCollapsedXiGetPtr
  END INTERFACE !CMISSBasisCollapsedXiGet
  
  !>Sets/changes the collapsed Xi flags for a basis.
  INTERFACE CMISSBasisCollapsedXiSet
    MODULE PROCEDURE CMISSBasisCollapsedXiSetNumber
    MODULE PROCEDURE CMISSBasisCollapsedXiSetPtr
  END INTERFACE !CMISSBasisCollapsedXiSet
  
  !>Finishes the creation of a new basis. \see OPENCMISS::CMISSBasisCreateStart
  INTERFACE CMISSBasisCreateFinish
    MODULE PROCEDURE CMISSBasisCreateFinishNumber
    MODULE PROCEDURE CMISSBasisCreateFinishPtr
  END INTERFACE !CMISSBasisCreateFinish
  
  !>Starts the creation of a new basis. \see OPENCMISS::CMISSBasisCreateFinish
  INTERFACE CMISSBasisCreateStart
    MODULE PROCEDURE CMISSBasisCreateFinishNumber
    MODULE PROCEDURE CMISSBasisCreateFinishPtr
  END INTERFACE !CMISSBasisCreateFinish
  
  !>Destroys a basis.
  INTERFACE CMISSBasisDestroy
    MODULE PROCEDURE CMISSBasisDestroyNumber
    MODULE PROCEDURE CMISSBasisDestroyPtr
  END INTERFACE !CMISSBasisDestroy

  !>Get the interpolation type in each Xi directions for a basis.
  INTERFACE CMISSBasisInterpolationXiGet
    MODULE PROCEDURE CMISSBasisInterpolationXiGetNumber
    MODULE PROCEDURE CMISSBasisInterpolaitonXiGetPtr
  END INTERFACE !CMISSBasisInterpolationXiGet
  
  !>Sets/changes the interpolation type in each Xi directions for a basis.
  INTERFACE CMISSBasisInterpolationXiSet
    MODULE PROCEDURE CMISSBasisInterpolationXiSetNumber
    MODULE PROCEDURE CMISSBasisInterpolaitonXiSetPtr
  END INTERFACE !CMISSBasisInterpolationXiSet
  
  !>Returns the number of local nodes in a basis.
  INTERFACE CMISSBasisNumberOfLocalNodesGet
    MODULE PROCEDURE CMISSBasisNumberOfLocalNodesGetNumber
    MODULE PROCEDURE CMISSBasisNumberOfLocalNodesGetPtr
  END INTERFACE !CMISSBasisNumberOfLocalNodesGet
  
  !>Returns the number of Xi directions in a basis.
  INTERFACE CMISSBasisNumberOfXiGet
    MODULE PROCEDURE CMISSBasisNumberOfXiGetNumber
    MODULE PROCEDURE CMISSBasisNumberOfXiGetPtr
  END INTERFACE !CMISSBasisNumberOfXiGet
  
  !>Sets/changes the number of Xi directions in a basis.
  INTERFACE CMISSBasisNumberOfXiSet
    MODULE PROCEDURE CMISSBasisNumberOfXiSetNumber
    MODULE PROCEDURE CMISSBasisNumberOfXiSetPtr
  END INTERFACE !CMISSBasisNumberOfXiSet
  
  !>Returns the number of Gauss points in each Xi direction on a basis quadrature.
  INTERFACE CMISSBasisQuadratureNumberOfGaussXiGet
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiGetNumber
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiGetPtr
  END INTERFACE !CMISSBasisQuadratureNumberOfGaussXiGet
  
  !>Sets/changes the number of Gauss points in each Xi direction on a basis quadrature.
  INTERFACE CMISSBasisQuadratureNumberOfGaussXiSet
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiSetNumber
    MODULE PROCEDURE CMISSBasisQuadratureNumberOfGaussXiSetPtr
  END INTERFACE !CMISSBasisQuadratureNumberOfGaussXiSet
  
  !>Returns the order of quadrature for a basis quadrature.
  INTERFACE CMISSBasisQuadratureOrderGet
    MODULE PROCEDURE CMISSBasisQuadratureOrderGetNumber
    MODULE PROCEDURE CMISSBasisQuadratureOrderGetPtr
  END INTERFACE !CMISSBasisQuadratureOrderGet
  
  !>Sets/changes the order of quadrature for a basis quadrature.
  INTERFACE CMISSBasisQuadratureOrderSet
    MODULE PROCEDURE CMISSBasisQuadratureOrderSetNumber
    MODULE PROCEDURE CMISSBasisQuadratureOrderSetPtr
  END INTERFACE !CMISSBasisQuadratureOrderSet
  
  !>Returns the quadrature type for a basis quadrature.
  INTERFACE CMISSBasisQuadratureTypeGet
    MODULE PROCEDURE CMISSBasisQuadratureTypeGetNumber
    MODULE PROCEDURE CMISSBasisQuadratureTypeGetPtr
  END INTERFACE !CMISSBasisQuadratureTypeGet
  
  !>Sets/changes the quadrature type for a basis quadrature.
  INTERFACE CMISSBasisQuadratureTypeSet
    MODULE PROCEDURE CMISSBasisQuadratureTypeSetNumber
    MODULE PROCEDURE CMISSBasisQuadratureTypeSetPtr
  END INTERFACE !CMISSBasisQuadratureTypeSet
  
  !>Returns the type of a basis.
  INTERFACE CMISSBasisTypeGet
    MODULE PROCEDURE CMISSBasisTypeGetNumber
    MODULE PROCEDURE CMISSBasisTypeGetPtr
  END INTERFACE !CMISSBasisTypeGet
  
  !>Sets/changes the type of a basis.
  INTERFACE CMISSBasisTypeSet
    MODULE PROCEDURE CMISSBasisTypeSetNumber
    MODULE PROCEDURE CMISSBasisTypeSetPtr
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
    MODULE PROCEDURE CMISSBoundaryConditionsDestroyPtr
  END INTERFACE !CMISSBoundaryConditionsDestroy

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant.
  INTERFACE CMISSBoundaryConditionsAddConstant
    MODULE PROCEDURE CMISSBoundaryConditionsAddConstantNumber
    MODULE PROCEDURE CMISSBoundaryConditionsAddConstantPtr
  END INTERFACE !CMISSBoundaryConditionsAddConstant

  !>Sets the value of the specified constant as a boundary condition on the specified constant.
  INTERFACE CMISSBoundaryConditionsSetConstant
    MODULE PROCEDURE CMISSBoundaryConditionsSetConstantNumber
    MODULE PROCEDURE CMISSBoundaryConditionsSetConstantPtr
  END INTERFACE !CMISSBoundaryConditionsSetConstant

  !>Adds to the value of the element constant and sets this as a boundary condition on the specified element.
  INTERFACE CMISSBoundaryConditionsAddElement
    MODULE PROCEDURE CMISSBoundaryConditionsAddElementNumber
    MODULE PROCEDURE CMISSBoundaryConditionsAddElementPtr
  END INTERFACE !CMISSBoundaryConditionsAddElement

  !>Sets the value of the specified element as a boundary condition on the specified element.
  INTERFACE CMISSBoundaryConditionsSetElement
    MODULE PROCEDURE CMISSBoundaryConditionsSetElementNumber
    MODULE PROCEDURE CMISSBoundaryConditionsSetElementPtr
  END INTERFACE !CMISSBoundaryConditionsSetElement

  !>Adds to the value of the node constant and sets this as a boundary condition on the specified node.
  INTERFACE CMISSBoundaryConditionsAddNode
    MODULE PROCEDURE CMISSBoundaryConditionsAddNodeNumber
    MODULE PROCEDURE CMISSBoundaryConditionsAddNodePtr
  END INTERFACE !CMISSBoundaryConditionsAddNode

  !>Sets the value of the specified node as a boundary condition on the specified node.
  INTERFACE CMISSBoundaryConditionsSetNode
    MODULE PROCEDURE CMISSBoundaryConditionsSetNodeNumber
    MODULE PROCEDURE CMISSBoundaryConditionsSetNodePtr
  END INTERFACE !CMISSBoundaryConditionsSetNode

  !>Gets the boundary conditions for an equations set.
  INTERFACE CMISSEquationsSetBoundaryConditionsGet
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsGetNumber
    MODULE PROCEDURE CMISSEquationsSetBoundaryConditionsGetPtr
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
    MODULE PROCEDURE CMISSControlLoopCurrentTimesGetPtr
  END INTERFACE !CMISSControlLoopCurrentTimesGet

  !>Destroy a control loop.
  INTERFACE CMISSControlLoopDestroy
    MODULE PROCEDURE CMISSControlLoopDestroyNumber0
    MODULE PROCEDURE CMISSControlLoopDestroyNumber1
    MODULE PROCEDURE CMISSControlLoopDestroyPtr
  END INTERFACE !CMISSControlLoopDestroy

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  INTERFACE CMISSControlLoopGet
    MODULE PROCEDURE CMISSControlLoopGetNumber00
    MODULE PROCEDURE CMISSControlLoopGetNumber10
    MODULE PROCEDURE CMISSControlLoopGetNumber01
    MODULE PROCEDURE CMISSControlLoopGetNumber11
    MODULE PROCEDURE CMISSControlLoopGetPtr0
    MODULE PROCEDURE CMISSControlLoopGetPtr1
  END INTERFACE !CMISSControlLoopGet

  !>Sets/changes the iteration parameters for a fixed control loop. \todo need a get metod
  INTERFACE CMISSControlLoopIterationsSet
    MODULE PROCEDURE CMISSControlLoopIterationsSetNumber0
    MODULE PROCEDURE CMISSControlLoopIterationsSetNumber1
    MODULE PROCEDURE CMISSControlLoopIterationsSetPtr
  END INTERFACE !CMISSControlLoopIterationsSet

  !>Sets/changes the maximum iterations for a while control loop. \todo need a get method
  INTERFACE CMISSControlLoopMaximumIterationsSet
    MODULE PROCEDURE CMISSControlLoopMaximumIterationsSetNumber0
    MODULE PROCEDURE CMISSControlLoopMaximumIterationsSetNumber1
    MODULE PROCEDURE CMISSControlLoopMaximumIterationsSetPtr
  END INTERFACE !CMISSControlLoopMaximumIterationsSet

  !>Returns the number of sub loops for a control loop.
  INTERFACE CMISSControlLoopNumberOfSubLoopsGet
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsGetNumber0
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsGetNumber1
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsGetPtr
  END INTERFACE !CMISSControlLoopNumberOfSubLoopsGet

  !>Sets/changes the number of sub loops for a control loop. \todo is this really a public method???
  INTERFACE CMISSControlLoopNumberOfSubLoopsSet
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsSetNumber0
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsSetNumber1
    MODULE PROCEDURE CMISSControlLoopNumberOfSubLoopsSetPtr
  END INTERFACE !CMISSControlLoopNumberOfSubLoopsGet

  !>Returns the time parameters for a time control loop.
  INTERFACE CMISSControlLoopTimesGet
    MODULE PROCEDURE CMISSControlLoopTimesGetNumber0
    MODULE PROCEDURE CMISSControlLoopTimesGetNumber1
    MODULE PROCEDURE CMISSControlLoopTimesGetPtr
  END INTERFACE !CMISSControlLoopTimesGet

  !>Sets/Changes the time parameters for a time control loop.
  INTERFACE CMISSControlLoopTimesSet
    MODULE PROCEDURE CMISSControlLoopTimesSetNumber0
    MODULE PROCEDURE CMISSControlLoopTimesSetNumber1
    MODULE PROCEDURE CMISSControlLoopTimesSetPtr
  END INTERFACE !CMISSControlLoopTimesSet

  !>Sets/Changes the loop type for a control loop. \todo Is this really a public method? \todo need a get method
  INTERFACE CMISSControlLoopTypeSet
    MODULE PROCEDURE CMISSControlLoopTypeSetNumber0
    MODULE PROCEDURE CMISSControlLoopTypeSetNumber1
    MODULE PROCEDURE CMISSControlLoopTypeSetPtr
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
    MODULE PROCEDURE CMISSCoordinateSystemCreateFinishPtr
  END INTERFACE !CMISSCoordinateSystemCreateFinish

  !>Starts the creation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemCreateFinish
  INTERFACE CMISSCoordinateSystemCreateStart
    MODULE PROCEDURE CMISSCoordinateSystemCreateStartNumber
    MODULE PROCEDURE CMISSCoordinateSystemCreateStartPtr
  END INTERFACE !CMISSCoordinateSystemCreateStart

  !>Destorys a coordinate system.
  INTERFACE CMISSCoordinateSystemDestroy
    MODULE PROCEDURE CMISSCoordinateSystemDestroyNumber
    MODULE PROCEDURE CMISSCoordinateSystemDestroyPtr
  END INTERFACE !CMISSCoordinateSystemDestroy

  !>Returns the coordinate system dimension. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemDimensionGet
    MODULE PROCEDURE CMISSCoordinateSystemDimensionGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemDimensionGetPtr
  END INTERFACE !CMISSCoordinateSystemDimensionGet

  !>Sets/changes the coordinate system dimension. \todo fix pointers
  INTERFACE CMISSCoordinateSystemDimensionSet
    MODULE PROCEDURE CMISSCoordinateSystemDimensionSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemDimensionSetPtr
  END INTERFACE !CMISSCoordinateSystemDimensionSet

  !>Returns the coordinate system focus. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemFocusGet
    MODULE PROCEDURE CMISSCoordinateSystemFocusGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemFocusGetPtr
  END INTERFACE !CMISSCoordinateSystemFocusGet
    
  !>Sets/changes the coordinate system focus. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemFocusSet
    MODULE PROCEDURE CMISSCoordinateSystemFocusSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemFocusSetPtr
  END INTERFACE !CMISSCoordinateSystemFocusSet

  !>Returns the coordinate system radial interpolation type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemRadialInterpolationGet
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationGetPtr
  END INTERFACE !CMISSCoordinateSystemRadialInterpolationGet
    
  !>Sets/changes the coordinate system radial interpolation type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemRadialInterpolationSet
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemRadialInterpolationSetPtr
  END INTERFACE !CMISSCoordinateSystemRadialInterpolationSet
    
  !>Returns the coordinate system type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemTypeGet
    MODULE PROCEDURE CMISSCoordinateSystemTypeGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemTypeGetPtr
  END INTERFACE !CMISSCoordinateSystemTypeGet
    
  !>Sets/changes the coordinate system type. \todo user number method \todo fix pointers
  INTERFACE CMISSCoordinateSystemTypeSet
    MODULE PROCEDURE CMISSCoordinateSystemTypeSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemTypeSetPtr
  END INTERFACE !CMISSCoordinateSystemTypeSet

  !>Returns the coordinate system orign. 
  INTERFACE CMISSCoordinateSystemOriginGet
    MODULE PROCEDURE CMISSCoordinateSystemOriginGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOriginGetPtr
  END INTERFACE !CMISSCoordinateSystemOriginGet

  !>Sets/changes the coordinate system orign. 
  INTERFACE CMISSCoordinateSystemOriginSet
    MODULE PROCEDURE CMISSCoordinateSystemOriginSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOriginSetPtr
  END INTERFACE !CMISSCoordinateSystemOriginSet

  !>Returns the coordinate system orientation. 
  INTERFACE CMISSCoordinateSystemOrientationGet
    MODULE PROCEDURE CMISSCoordinateSystemOrientationGetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOrientationGetPtr
  END INTERFACE !CMISSCoordinateSystemOrientationGet

  !>Sets/changes the coordinate system orientation. 
  INTERFACE CMISSCoordinateSystemOrientationSet
    MODULE PROCEDURE CMISSCoordinateSystemOrientationSetNumber
    MODULE PROCEDURE CMISSCoordinateSystemOrientationSetPtr
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
  INTEGER(INTG), PARAMETER :: CMISSFieldDelUDelNVariableType = FIELD_DELUDELN_VARIABLE_TYPE !<Normal derivative variable type i.e., du/dn \see OPENCMISS_FieldVariableTypes,OEPNCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDelUDelTVariableType = FIELD_DELUDELT_VARIABLE_TYPE !<First time derivative variable type i.e., du/dt \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDel2UDelT2VariableType = FIELD_DEL2UDELT2_VARIABLE_TYPE !<Second type derivative variable type i.e., d^2u/dt^2 \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldVVariableType = FIELD_V_VARIABLE_TYPE !<Second standard variable type i.e., v \see OPENCMISS_FieldVariableTypes,OPENCMISS
  INTEGER(INTG), PARAMETER :: CMISSFieldDelVDelNVariableType = FIELD_DELVDELN_VARIABLE_TYPE !<Second normal variable type i.e., dv/dn \see OPENCMISS_FieldVariableTypes,OPENCMISS
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSFieldDependentType,CMISSFieldIndependentType

  PUBLIC CMISSFieldScalarDimensionType,CMISSFieldVectorDimensionType,CMISSFieldTensorDimensionType

  PUBLIC CMISSFieldGeometricType,CMISSFieldFibreType,CMISSFieldGeneralType,CMISSFieldMaterialType

  PUBLIC CMISSFieldConstantInterpolation,CMISSFieldElementBasedInterpolation,CMISSFieldNodeBasedInterpolation, &
    & CMISSFieldGridPointBasedInterpolation,CMISSFieldGaussPointBasedInterpolation

  PUBLIC CMISSFieldUVariableType,CMISSFieldDelUDelNVariableType,CMISSFieldDelUDelTVariableType,CMISSFieldDel2UDelT2VariableType, &
    & CMISSFieldVVariableType,CMISSFieldDelVDelNVariableType

  
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
  INTEGER(INTG), PARAMETER :: CMISSProblemReactionDiffusionEquationType = PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE=7 !<Reaction-Diffusion equation problem type \see OPENCMISS_ProblemTypes,OPENCMISS 
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
    CALL CMISS_Initialise(WORLD_COORDINATE_SYSTEM,WORLD_REGION,Err,ERROR,*999)
    WorldCoordianteUserNumber=WORLD_COORDINATE_SYSTEM%USER_NUMBER
    WorldRegionUserNumber=WORLD_REGION%USER_NUMBER

    RETURN
999 CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSInitialiseNumber

  !
  !================================================================================================================================
  !
  
  !>Initialises CMISS returning a pointer to the world coordinate system and region.
  SUBROUTINE CMISSInitialisePtr(WorldCoordinateSystem,WorldRegion,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(OUT), POINTER :: WorldCoordinateSystem !<On return, a pointer to the world coordinate system.
    TYPE(REGION_TYPE), INTENT(OUT), POINTER :: WorldRegion !<On return, a pointer to the world region.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: WORLD_COORDINATE_SYSTEM
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION

    CALL CMISS_Initialise(WorldCoordinateSystem,WorldRegion,Err,ERROR,*999)

    RETURN
999 CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSInitialisePtr


!!==================================================================================================================================
!!
!! ANALYTIC_ANALYSIS_ROUTINES
!!
!!==================================================================================================================================

  !>Output the analytic error analysis for a field specified by a user number compared to the analytic values parameter set.
  SUBROUTINE CMISSAnalyticAnalytisOutputNumber(UserNumber,FileName,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the field to calculate the analytic error analysis for.
    CHARACTER(LEN=*) :: FileName !<If not empty, the filename to output the analytic analysis to. If empty, the analysis will be output to the standard output.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(FIELD_TYPE), POINTER :: FIELD

    CALL ENTERS("CMISSAnalyticAnalytisOutputNumber",Err,ERROR,*999)
    
    NULLIFY(FIELD)
    CALL FIELD_USER_NUMBER_FIND(UserNumber,FIELD,Err,ERROR,*999)
    CALL ANALYTIC_ANALYSIS_OUTPUT(Field,FileName,Err,ERROR,*999)
    
    CALL EXITS("CMISSAnalyticAnalytisOutputNumber")
    RETURN
999 CALL ERRORS("CMISSAnalyticAnalytisOutputNumber",Err,ERROR)
    CALL EXITS("CMISSAnalyticAnalytisOutputNumber")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSAnalyticAnalysisOutputNumber

  !
  !================================================================================================================================
  !  

  !>Output the analytic error analysis for a field identified by a pointer compared to the analytic values parameter set.
  SUBROUTINE CMISSAnalyticAnalytisOutputPtr(Field,FileName,Err)
  
    !Argument variables
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: Field !<A pointer to the dependent field to calculate the analytic error analysis for.
    CHARACTER(LEN=*) :: FileName !<If not empty, the filename to output the analytic analysis to. If empty, the analysis will be output to the standard output.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSAnalyticAnalytisOutputPtr",Err,ERROR,*999)
    
    CALL ANALYTIC_ANALYSIS_OUTPUT(Field,FileName,Err,ERROR,*999)

    CALL EXITS("CMISSAnalyticAnalytisOutputPtr")
    RETURN
999 CALL ERRORS("CMISSAnalyticAnalytisOutputPtr",Err,ERROR)
    CALL EXITS("CMISSAnalyticAnalytisOutputPtr")    
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSAnalyticAnalysisOutputPtr

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

    CALL ENTERS("CMISSBasisCollapsedXiGetNumber",Err,ERROR,*999)

    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_COLLAPSED_XI_GET(BASIS,CollapsedXi,Err,ERROR,*999)

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
  
  !>Returns the collapsed Xi flags of a basis identified by a pointer.
  SUBROUTINE CMISSBasisCollapsedXiGetPtr(Basis,CollapsedXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the collapsed Xi flags for.
    INTEGER(INTG), INTENT(OUT) :: CollapsedXi(:) !<CollapsedXi(ni). On return, the collapsed Xi parameter for the ni'th Xi direction. \see OPENCMISS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisCollapsedXiGetPtr",Err,ERROR,*999)
    
    CALL BASIS_COLLAPSED_XI_GET(Basis,CollapsedXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCollapsedXiGetPtr"
    RETURN
999 CALL ERRORS("CMISSBasisCollapsedXiGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisCollapsedXiGetPtr"
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCollapsedXiGetPtr

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

    CALL ENTERS("CMISSBasisCollapsedXiSetNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_COLLAPSED_XI_SET(BASIS,CollapsedXi,Err,ERROR,*999)

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
  
  !>Sets/changes the collapsed Xi flags of a basis identified by a pointer.
  SUBROUTINE CMISSBasisCollapsedXiGetPtr(Basis,CollapsedXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to set the collapsed Xi flags for.
    INTEGER(INTG), INTENT(IN) :: CollapsedXi(:) !<CollapsedXi(ni). The collapsed Xi parameter for the ni'th Xi direction to set. \see OPENCMISS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisCollapsedXiGetPtr",Err,ERROR,*999)

    CALL BASIS_COLLAPSED_XI_SET(Basis,CollapsedXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCollapsedXiGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisCollapsedXiGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisCollapsedXiGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCollapsedXiSetPtr

  !
  !================================================================================================================================
  !
  
   !>Finishes the creation of a new basis identified by a user number.
  SUBROUTINE CMISSBasisCreateFinishNumber(UserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to finish the creation of.
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS

    CALL ENTERS("CMISSBasisCreateFinishNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,Err,ERROR,*999)
    CALL BASIS_CREATE_FINISH(Basis,Err,ERROR,*999)

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
  
  !>Finishes the creation of a new basis identified by a pointer.
  SUBROUTINE CMISSBasisCreateFinishPtr(Basis,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code     
    !Local variables

    CALL ENTERS("CMISSBasisCreateFinishPtr",Err,ERROR,*999)

    CALL BASIS_CREATE_FINISH(Basis,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCreateFinishPtr")
    RETURN
999 CALL ERRORS("CMISSBasisCreateFinishPtr",Err,ERROR)
    CALL EXITS("CMISSBasisCreateFinishPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCreateFinishPtr
  
  !
  !================================================================================================================================
  !
  
  !>Starts the creation of a new basis.
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
 
  !>Starts the creation of a new basis returning a pointer.
  SUBROUTINE CMISSBasisCreateStartPtr(UserNumber,Basis,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to start the creation of.
    TYPE(BASIS_TYPE), POINTER :: Basis !<On exit, a pointer to the basis to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisCreateStartPtr",Err,ERROR,*999)
    
    CALL BASIS_CREATE_START(UserNumber,Basis,Err,ERROR,*999)

    CALL EXITS("CMISSBasisCreateStartPtr")
    RETURN
999 CALL ERRORS("CMISSBasisCreateStartPtr",Err,ERROR)
    CALL EXITS("CMISSBasisCreateStartPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisCreateStartPtr
  
  !
  !================================================================================================================================
  !
  
  !>Destroys a basis identified by its basis user number.
  SUBROUTINE CMISSBasisDestroyNumber(UserNumber,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: UserNumber !<The user number of the basis to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(BASIS_TYPE), POIINTER :: BASIS

    CALL ENTERS("CMISSBasisDestroyNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_DESTROY(UserNumber,Err,ERROR,*999)

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
  
  !>Destroys a basis identified by a pointer.
  SUBROUTINE CMISSBasisDestroyPtr(Basis,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to destroy.
    !Local variables

    CALL ENTERS("CMISSBasisDestroyPtr",Err,ERROR,*999)
    
    CALL BASIS_DESTROY(Basis,Err,ERROR,*999)

    CALL EXITS("CMISSBasisDestroyPtr")
    RETURN
999 CALL ERRORS("CMISSBasisDestroyPtr",Err,ERROR)RETURN
    CALL EXITS("CMISSBasisDestroyPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisDestroyPtr
  
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
    
    CALL ENTERS("CMISSBasisInterpolationXiGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_GET(BASIS,Err,ERROR,*999)

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
  
  !>Get the interpolation type in each xi directions for a basis indentified by a pointer. \todo User number form
  SUBROUTINE CMISSBasisInterpolationXiGetPtr(Basis,InterpolationXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the interpolation xi for.
    INTEGER(INTG), INTENT(OUT) :: InterpolationXi(:) !<On return, the interpolation xi parameters for each Xi direction \see OPENCMISS_InterpolationSpecifications.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisInterpolationXiGetPtr",Err,ERROR,*999)
    
    CALL BASIS_INTERPOLATION_XI_GET(Basis,Err,ERROR,*999)

    CALL EXITS("CMISSBasisInterpolationXiGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisInterpolationXiGetPtr",Err,ERROR)
    CALL CMISS_HANDL_ERROR(Err,ERROR)
    CALL EXITS("CMISSBasisInterpolationXiGetPtr")
    RETURN
    
  END SUBROUTINE CMISSBasisInterpolationXiGetPtr
  
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

    CALL ENTERS("CMISSBasisInterpolationXiSetNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_SET(BASIS,Err,ERROR,*999)

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
  
  !>Sets/changes the interpolation type in each xi directions for a basis indentified by a pointer.
  SUBROUTINE CMISSBasisInterpolationXiGetPtr(Basis,InterpolationXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the interpolation xi for.
    INTEGER(INTG), INTENT(IN) :: InterpolationXi(:) !<The interpolation xi parameters for each Xi direction \see OPENCMISS_InterpolationSpecifications.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisInterpolationXiGetPtr",Err,ERROR,*999)
    
    CALL BASIS_INTERPOLATION_XI_SET(Basis,Err,ERROR,*999)

    CALL EXITS("CMISSBasisInterpolationXiGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisInterpolationXiGetPtr",Err,ERORR)
    CALL EXITS("CMISSBasisInterpolationXiGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisInterpolationXiSetPtr
  
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

    CALL ENTERS("CMISSBasisNumberOfLocalNodesGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_NUMBER_OF_LOCAL_NODES_GET(BASIS,NumberOfLocalNodes,Err,ERROR,*999)

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
  
  !>Returns the number of local nodes in a basis identified by a pointer.
  SUBROUTINE CMISSBasisNumberOfLocalNodesGetPtr(Basis,NumberOfLocalNodes,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the number of local nodes for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfLocalNodes !<On return, the number of local nodes in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisNumberOfLocalNodesGetPtr",Err,ERROR,*999)

    CALL BASIS_NUMBER_OF_LOCAL_NODES_GET(Basis,NumberOfLocalNodes,Err,ERROR,*999)

    CALL EXITS("CMISSBasisNumberOfLocalNodesGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfLocalNodesGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfLocalNodesGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfLocalNodesGetPtr
  
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

    CALL ENTERS("CMISSBasisNumberOfXiGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_NUMBER_OF_XI_GET(BASIS,NumberOfXi,Err,ERROR,*999)

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
  
  !>Returns the number of Xi directions in a basis identified by a pointer.
  SUBROUTINE CMISSBasisNumberOfXiGetPtr(Basis,NumberOfXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the number of xi directions for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfXi !<On return, the number of xi directions in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisNumberOfXiGetPtr",Err,ERROR,*999)
    
    CALL BASIS_NUMBER_OF_XI_GET(BASIS,NumberOfXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisNumberOfXiGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfXiGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfXiGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfXiGetPtr
  
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

    CALL ENTERS("CMISSBasisNumberOfXiSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_NUMBER_OF_XI_SET(BASIS,NumberOfXi,Err,ERROR,*999)

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
  
  !>Sets/changes the number of Xi directions in a basis identified by a pointer.
  SUBROUTINE CMISSBasisNumberOfXiSetPtr(Basis,NumberOfXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to set the number of xi directions for.
    INTEGER(INTG), INTENT(IN) :: NumberOfXi !<The number of xi directions in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisNumberOfXiSetPtr",Err,ERROR,*999)
    
    CALL BASIS_NUMBER_OF_XI_SET(BASIS,NumberOfXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisNumberOfXiSetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisNumberOfXiSetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisNumberOfXiSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisNumberOfXiSetPtr
  
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

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET(BASIS,NumberOfGaussXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetNumber",Err,ERROR,*999)
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiGetNumber",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetNumber",Err,ERROR,*999)
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetNumber
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number Gauss points in each Xi directions for a basis quadrature identified by a pointer.
  SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetPtr(Basis,NumberOfGaussXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the number of Gauss Xi for.
    INTEGER(INTG), INTENT(OUT) :: NumberOGaussfXi(:) !<On return, the number of Gauss points in each Xi directions in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiGetPtr",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET(BASIS,NumberOfGaussXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiGetPtr
  
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

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(BASIS,NumberOfGaussXi,Err,ERROR,*999)

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
  
  !>Sets the number Gauss points in each Xi directions for a basis quadrature identified by a pointer.
  SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiSetPtr(Basis,NumberOfGaussXi,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the number of Gauss Xi for.
    INTEGER(INTG), INTENT(IN) :: NumberOGaussfXi(:) !<The number of Gauss points in each Xi directions in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureNumberOfGaussXiSetPtr",Err,ERROR,*999)

    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(BASIS,NumberOfGaussXi,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiSetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureNumberOfGaussXiSetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureNumberOfGaussXiSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureNumberOfGaussXiSetPtr
  
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

    CALL ENTERS("CMISSBasisQuadratureOrderGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_ORDER_GET(BASIS,QuadratureOrder,Err,ERROR,*999)

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
  
  !>Returns the the order of quadrature for a basis quadrature identified by a pointer.
  SUBROUTINE CMISSBasisQuadratureOrderGetPtr(Basis,QuadratureOrder,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), INTENT(IN), POINTER :: Basis !<A pointer to the basis to get the quadrature order for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureOrder !<On return, the order of quadrature in the specified basis.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureOrderGetPtr",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_ORDER_GET(BASIS,QuadratureOrder,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureOrderGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureOrderGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureOrderGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureOrderGetPtr

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

    CALL ENTERS("CMISSBasisQuadratureOrderSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_ORDER_SET(BASIS,QuadratureOrder,Err,ERROR,*999)

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
  
  !>Sets/changes the the order of quadrature for a basis quadrature identified by a pointer.
  SUBROUTINE CMISSBasisQuadratureOrderSetPtr(Basis,QuadratureOrder,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to set the quadrature order for.
    INTEGER(INTG), INTENT(IN) :: QuadratureOrder !<The order of quadrature in the specified basis to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureOrderSetPtr",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_ORDER_SET(BASIS,QuadratureOrder,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureOrderSetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureOrderSetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureOrderSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureOrderSetPtr
  
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

    CALL ENTERS("CMISSBasisQuadratureTypeGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_TYPE_GET(BASIS,QuadratureType,Err,ERROR,*999)

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
  
  !>Returns the the type of quadrature for a basis quadrature identified by a pointer.
  SUBROUTINE CMISSBasisQuadratureTypeGetPtr(Basis,QuadratureType,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the quadrature order for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureType !<On return, the type of quadrature in the specified basis. \see OPENCMISS_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureTypeGetPtr",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_TYPE_GET(BASIS,QuadratureType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureTypeGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureTypeGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureTypeGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureTypeGetPtr

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

    CALL ENTERS("CMISSBasisQuadratureTypeSetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_TYPE_SET(BASIS,QuadratureType,Err,ERROR,*999)

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
  
  !>Sets/changes the the type of quadrature for a basis quadrature identified by a pointer.
  SUBROUTINE CMISSBasisQuadratureTypeSetPtr(Basis,QuadratureType,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the quadrature type for.
    INTEGER(INTG), INTENT(OUT) :: QuadratureType !<The type of quadrature in the specified basis to set. \see OPENCMISS_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisQuadratureTypeSetPtr",Err,ERROR,*999)
    
    CALL BASIS_QUADRATURE_TYPE_SET(Basis,QuadratureType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisQuadratureTypeSetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisQuadratureTypeSetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisQuadratureTypeSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisQuadratureTypeSetPtr

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

    CALL ENTERS("CMISSBasisTypeGetNumber",Err,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_TYPE_GET(BASIS,BasisType,Err,ERROR,*999)

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
  
  !>Returns the type of a basis identified by a pointer.
  SUBROUTINE CMISSBasisTypeGetPtr(Basis,BasisType,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to get the type for.
    INTEGER(INTG), INTENT(OUT) :: BasisType !<On return, the type of the specified basis. \see OPENCMISS_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisTypeGetPtr",Err,ERROR,*999)
    
    CALL BASIS_TYPE_GET(Basis,BasisType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisTypeGetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisTypeGetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisTypeGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeGetPtr

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

    CALL ENTERS("CMISSBasisTypeSetNumber",ERR,ERROR,*999)
    
    NULLIFY(BASIS)
    CALL BASIS_USER_NUMBER_FIND(UserNumber,BASIS,ERR,ERROR,*999)
    CALL BASIS_TYPE_SET(BASIS,BasisType,Err,ERROR,*999)

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
  
  !>Sets/changes the type of a basis identified by a pointer.
  SUBROUTINE CMISSBasisTypeSetPtr(Basis,BasisType,Err)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: Basis !<A pointer to the basis to set the type for.
    INTEGER(INTG), INTENT(IN) :: BasisType !<The type of the specified basis to set. \see OPENCMISS_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBasisTypeSetPtr",Err,ERROR,*999)
    
    CALL BASIS_TYPE_SET(Basis,BasisType,Err,ERROR,*999)

    CALL EXITS("CMISSBasisTypeSetPtr")
    RETURN
999 CALL ERRORS("CMISSBasisTypeSetPtr",Err,ERROR)
    CALL EXITS("CMISSBasisTypeSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBasisTypeSetPtr


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

    CALL ENTERS("CMISSBoundaryConditionsDestroyNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_DESTROY(BOUNDARY_CONDITIONS,Err,ERROR,*999)

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
    
  !>Destroys boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsDestroyPtr(BoundaryConditions,Err)
  
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsDestroyPtr",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_DESTROY(BoundaryConditions,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsDestroyPtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsDestroyPtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsDestroyPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsDestroyPtr

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

    CALL ENTERS("CMISSBoundaryConditionsAddConstantNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_ADD_CONSTANT(BOUNDARY_CONDITIONS,VariableType,ComponentNumber,Condition,Value,Err,ERROR,*999)

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
  
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant for boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsAddConstantPtr(BoundaryConditions,VariableType,ComponentNumber,Condition,Value,Err)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to add the constant to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSBoundaryConditionsAddConstantPtr",Err,ERROR,*999)
  
    CALL BOUNDARY_CONDITIONS_ADD_CONSTANT(BoundaryConditions,VariableType,ComponentNumber,Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsAddConstantPtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddConstantPtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddConstantPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddConstantPtr
 
  !
  !================================================================================================================================
  !
  

  !>Sets the value of the specified constant as a boundary condition on the specified constant for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsSetConstantNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,ComponentNumber, &
    &  Condition,Value,Err)
  
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

    CALL ENTERS("CMISSBoundaryConditionsSetConstantNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_SET_CONSTANT(BOUNDARY_CONDITIONS,VariableType,ComponentNumber,Condition,Value,Err,ERROR,*999)

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
  
  !>Sets the value of the specified constant and sets this as a boundary condition on the specified constant for boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsSetConstantPtr(BoundaryConditions,VariableType,ComponentNumber,Condition,Value,Err)
  
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to set the constant to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSBoundaryConditionsSetConstantPtr",ERR,ERROR,*999)
  
    CALL BOUNDARY_CONDITIONS_SET_CONSTANT(BoundaryConditions,VariableType,ComponentNumber,Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsSetConstantPtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetConstantPtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetConstantPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetConstantPtr

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

    CALL ENTERS("CMISSBoundaryConditionsAddElementNumber",Err,ERROR,*999)

    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_ADD_ELEMENT(BOUNDARY_CONDITIONS,VariableType,ElementUserNumber,ComponentNumber,Condition,Value, &
      & Err,ERROR,*999)

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
  
  !>Adds to the value of the specified element and sets this as a boundary condition on the specified element for boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsAddElementPtr(BoundaryConditions,VariableType,ElementUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to add the element to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ElementUserNumber !<The user number of the element to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsAddElementPtr",Err,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_ADD_ELEMENT(BoundaryConditions,VariableType,ElementUserNumber,ComponentNumber,Condition,Value, &
      & Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsAddElementPtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddElementPtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddElementPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddElementPtr
 
  !
  !================================================================================================================================
  !

  !>Sets the value of the specified element as a boundary condition on the specified element for boundary conditions identified by a user number.
  SUBROUTINE CMISSBoundaryConditionsSetElementNumber(RegionUserNumber,EquationsSetUserNumber,VariableType,ElementUserNumber, &
    & ComponentNumber,ElementUserNumber,Condition,Value,Err)
  
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

    CALL ENTERS("CMISSBoundaryConditionsSetElementNumber",Err,ERROR,*999)

    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_SET_ELEMENT(BOUNDARY_CONDITIONS,VariableType,ElementUserNumber,ComponentNumber,Condition,Value, &
      & Err,ERROR,*999)

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
  
  !>Sets the value of the specified element and sets this as a boundary condition on the specified elements for boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsSetElementPtr(BoundaryConditions,VariableType,ElementUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to set the element to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: ElementUserNumber !<The user number of the element to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsSetElementPtr",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_SET_ELEMENT(BoundaryConditions,VariableType,ElementUserNumber,ComponentNumber,Condition,Value, &
      & Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsSetElementPtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetElementPtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetElementPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetElementPtr

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

    CALL ENTERS("CMISSBoundaryConditionsAddNodeNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_ADD_NODE(BOUNDARY_CONDITIONS,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
      & Condition,Value,Err,ERROR,*999)

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
  
  !>Adds to the value of the specified node and sets this as a boundary condition on the specified node for boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsAddNodePtr(BoundaryConditions,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to add the node to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to add the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: DerivativeNumber !<The user number of the node derivative to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: NodeUserNumber !<The user number of the node to add the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to add.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSBoundaryConditionsAddNodePtr",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_ADD_NODE(BoundaryConditions,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
      & Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsAddNodePtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsAddNodePtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsAddNodePtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsAddNodePtr
 
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

    CALL ENTERS("CMISSBoundaryConditionsSetNodeNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(BOUNDARY_CONDITIONS)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,VariableType,DerivativeNumber,NodeUserNumber,ComponentNumber, &
      & Condition,Value,Err,ERROR,*999)

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
  
  !>Sets the value of the specified node and sets this as a boundary condition on the specified node for boundary conditions identified by a pointer.
  SUBROUTINE CMISSBoundaryConditionsSetElementPtr(BoundaryConditions,VariableType,DerivativeNumber,NodeUserNumberComponentNumber, &
    & Condition,Value,Err)
  
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BoundaryConditions !<A pointer to the boundary conditions to set the node to.
    INTEGER(INTG), INTENT(IN) :: VariableType !<The variable type of the dependent field to set the boundary condition at. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: DerivativeNumber !<The user number of the node derivative to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: NodeUserNumber !<The user number of the node to set the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: ComponentNumber !<The component number of the dependent field to set the boundary condition at.
    INTEGER(INTG), INTENT(IN) :: Condition !<The boundary condition type to set \see OPENCMISS_BoundaryConditions,OPENCMISS
    REAL(DP), INTENT(IN) :: Value !<The value of the boundary condition to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
  
    CALL ENTERS("CMISSBoundaryConditionsSetElementPtr",Err,ERROR,*999)
    
    CALL BOUNDARY_CONDITIONS_SET_NODE(BoundaryConditions,VariableType,DerivativeNumber,NodeUserNumberComponentNumber, &
      & Condition,Value,Err,ERROR,*999)

    CALL EXITS("CMISSBoundaryConditionsSetElementPtr")
    RETURN
999 CALL ERRORS("CMISSBoundaryConditionsSetElementPtr",Err,ERROR)
    CALL EXITS("CMISSBoundaryConditionsSetElementPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSBoundaryConditionsSetNodePtr

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for an equations set identified by a user number. 
  SUBROUTINE CMISSEquationsSetBoundaryConditionsGetNumber(RegionUserNumber,EquationsSetUserNumber,BoundaryConditions,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: RegionUserNumber !<The user number of the region containing the equations set to get the boundary conditions for.
    INTEGER(INTG), INTENT(IN) :: EquationsSetUserNumber !<The user number of the equations set to get the boundary conditions for.
    TYPE(BOUNDARY_CONDITIONS_TYPE), INTENT(OUT), POINTER :: BoundaryConditions, !<On return, a pointer to the boundary conditions for the specified equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(REGION_TYPE), POINTER :: REGION

    CALL ENTERS("CMISSEquationsSetBoundaryConditionsGetNumber",Err,ERROR,*999)
    
    NULLIFY(REGION)
    NULLIFY(EQUATIONS_SET)
    CALL REGION_USER_NUMBER_FIND(RegionUserNumber,REGION,Err,ERROR,*999)
    CALL EQUATIONS_SET_USER_NUMBER_FIND(EquationsSetUserNumber,REGION,EQUATIONS_SET,Err,ERROR,*999)
    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,Err,ERROR,*999)

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
  SUBROUTINE CMISSEquationsSetBoundaryConditionsGetPtr(EquationsSet,BoundaryConditions,Err)
  
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), INTENT(IN), POINTER :: EquationsSet !<A pointer to the equations set to get the boundary conditions for.
    TYPE(BOUNDARY_CONDITIONS_TYPE), INTENT(OUT), POINTER :: BoundaryConditions, !<On return, a pointer to the boundary conditions for the specified equations set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    
    CALL ENTERS("CMISSEquationsSetBoundaryConditionsGetPtr",ERR,ERROR,*999)

    CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EquationsSet,BoundaryConditions,Err,ERROR,*999)

    CALL EXITS("CMISSEquationsSetBoundaryConditionsGetPtr")
    RETURN
999 CALL ERRORS("CMISSEquationsSetBoundaryConditionsGetPtr",Err,ERROR)
    CALL EXITS("CMISSEquationsSetBoundaryConditionsGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSEquationsSetBoundaryConditionsGetPtr

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

    CALL ENTERS("CMISSControlLoopCurrentTimesGetNumber0",Err,ERROR,*999)
 
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CurrentTime,TimeIncrement,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopCurrentTimesGetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CurrentTime,TimeIncrement,Err,ERROR,*999)

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
  
  !>Gets the current time parameters for a time control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopCurrentTimesGetPtr(ControlLoop,CurrentTime,TimeIncrement,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to get the current times for.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time of the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the current time increment of the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopCurrentTimesGetPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_CURRENT_TIMES_GET(ControlLoop,CurrentTime,TimeIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopCurrentTimesGetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopCurrentTimesGetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopCurrentTimesGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopCurrentTimesGetPtr

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

    CALL ENTERS("CMISSControlLoopDestroyNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_DESTROY(CONTROL_LOOP,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopDestroyNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_DESTROY(CONTROL_LOOP,Err,ERROR,*999)

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
  
  !>Destroys a control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopDestroyPtr(ControlLoop,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(INOUT), POINTER :: ControlLoop !<A pointer to the control loop to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopDestroyPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_DESTROY(ControlLoop,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopDestroyPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopDestroyPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopDestroyPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopDestroyPtr
 
  !
  !================================================================================================================================
  !  
  
  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root identified by user numbers.
  SUBROUTINE CMISSControlLoopGetNumber00(ProblemUserNumber,ControlLoopRootIdentifer,ControlLoopIdentifier,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifier !<The root control loop identifier.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    TYPE(CONTROL_LOOP_TYPE), INTENT(OUT), POINTER :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("CMISSControlLoopGetNumber00",Err,ERROR,*999)
    
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifier,ROOT_CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifier,ControlLoop,Err,ERROR,*999)

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
  SUBROUTINE CMISSControlLoopGetNumber10(ProblemUserNumber,ControlLoopRootIdentifers,ControlLoopIdentifier,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifiers(:) !<The root control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    TYPE(CONTROL_LOOP_TYPE), INTENT(OUT), POINTER :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("CMISSControlLoopGetNumber10",Err,ERROR,*999)
 
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifiers,ROOT_CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifier,ControlLoop,Err,ERROR,*999)

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
  SUBROUTINE CMISSControlLoopGetNumber01(ProblemUserNumber,ControlLoopRootIdentifer,ControlLoopIdentifiers,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifier !<The root control loop identifier.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    TYPE(CONTROL_LOOP_TYPE), INTENT(OUT), POINTER :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("CMISSControlLoopGetNumber01",Err,ERROR,*999)
    
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifier,ROOT_CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifiers,ControlLoop,Err,ERROR,*999)

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
  SUBROUTINE CMISSControlLoopGetNumber11(ProblemUserNumber,ControlLoopRootIdentifers,ControlLoopIdentifiers,ControlLoop,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: ControlLoopRootIdentifiers(:) !<The root control loop identifiers.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    TYPE(CONTROL_LOOP_TYPE), INTENT(OUT), POINTER :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ROOT_CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("CMISSControlLoopGetNumber11",Err,ERROR,*999)
    
    NULLIFY(ROOT_CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopRootIdentifiers,ROOT_CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_GET(ROOT_CONTROL_LOOP,ControlLoopIdentifiers,ControlLoop,Err,ERROR,*999)

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
  
  !>Destroys a control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopGetPtr0(ControlLoopRoot,ControlLoopIdentifier,ControlLoop,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoopRoot !<A pointer to the root control loop.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifier !<The control loop identifier.
    TYPE(CONTROL_LOOP_TYPE), INTENT(OUT), POINTER :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopGetPtr0",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_GET(ControlLoopRoot,ControlLoopIdentifier,ControlLoop,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopGetPtr0")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetPtr0",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetPtr0")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetPtr0
  
  !
  !================================================================================================================================
  !  
  
  !>Destroys a control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopGetPtr1(ControlLoopRoot,ControlLoopIdentifiers,ControlLoop,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoopRoot !<A pointer to the root control loop.
    INTEGER(INTG), INTENT(IN) :: ControlLoopIdentifiers(:) !<The control loop identifiers.
    TYPE(CONTROL_LOOP_TYPE), INTENT(OUT), POINTER :: ControlLoop !<On return, the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopGetPtr1",Err,ERROR,*999)

    CALL CONTROL_LOOP_GET(ControlLoopRoot,ControlLoopIdentifiers,ControlLoop,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopGetPtr1")
    RETURN
999 CALL ERRORS("CMISSControlLoopGetPtr1",Err,ERROR)
    CALL EXITS("CMISSControlLoopGetPtr1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopGetPtr1

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

    CALL ENTERS("CMISSControlLoopIterationsSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_ITERATIONS_SET(CONTROL_LOOP,StartIteration,StopIteration,IterationIncrement,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopIterationsSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_ITERATIONS_SET(CONTROL_LOOP,StartIteration,StopIteration,IterationIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopIterationsSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopIterationsSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopIterationsSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopIterationsSetNumber1

  !================================================================================================================================
  !  
  
  !>Sets/changes the iteration parameters for a fixed control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopIterationsSetPtr(ControlLoop,StartIteration,StopIteration,IterationIncrement,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to set the iteration parameters for.
    INTEGER(INTG), INTENT(IN) :: StartIteration !<The start iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: StopIteration !<The stop iteration of the fixed control loop to set.
    INTEGER(INTG), INTENT(IN) :: IterationIncrement !<The iteration increment of the fixed control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopIterationsSetPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_ITERATIONS_SET(ControlLoop,StartIteration,StopIteration,IterationIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopIterationsSetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopIterationsSetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopIterationsSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopIterationsSetPtr

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

    CALL ENTERS("CMISSControlLoopMaximumIterationsSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MaximumIterations,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopMaximumIterationsSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MaximumIterations,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopMaximumIterationsSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopMaximumIterationsSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopMaximumIterationsSetNumber1")
    CALLL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopMaximumIterationsSetNumber1

  !
  !================================================================================================================================
  !  
  
  !>Sets/changes the maximum iterations for a while control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopMaximumIterationsSetPtr(ControlLoop,MaximumIterations,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to set the maximum iterations for.
    INTEGER(INTG), INTENT(IN) :: MaximumIterations !<The maximum iterations of the while control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopMaximumIterationsSetPtr",Err,ERROR)
    
    CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MaximumIterations,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopMaximumIterationsSetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopMaximumIterationsSetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopMaximumIterationsSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopMaximumIterationsSetPtr1

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

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsGetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsGetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)

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
  
  !>Returns the number of sub-control loops for a control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetPtr(ControlLoop,NumberOfSubLoops,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to get the number of sub loops for.
    INTEGER(INTG), INTENT(OUT) :: NumberOfSubLoops !<On return, the number of sub loops for the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsGetPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(ControlLoop,NumberOfSubLoops,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsGetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsGetPtr

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

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsSetNumber",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetNumber0")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsSetNumber0")
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

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,NumberOfSubLoops,Err,ERROR,*999)

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
  
  !>Sets/changes the number of sub-control loops for a control loop identified by a pointer. \todo is this really public???
  SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetPtr(ControlLoop,NumberOfSubLoops,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to set the number of sub loops for.
    INTEGER(INTG), INTENT(IN) :: NumberOfSubLoops !<The number of sub loops for the specified control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopNumberOfSubLoopsSetPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(ControlLoop,NumberOfSubLoops,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopNumberOfSubLoopsSetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopNumberOfSubLoopsSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopNumberOfSubLoopsSetPtr

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

    CALL ENTERS("CMISSControlLoopTimesGetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,CurrentTime,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopTimesGetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,CurrentTime,Err,ERROR,*999)

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
 
  !>Returns the time parameters for a time control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopTimesGetPtr(ControlLoop,StartTime,StopTime,TimeIncrement,CurrentTime,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to get the times for.
    REAL(DP), INTENT(OUT) :: StartTime !<On return, the start time for the time control loop.
    REAL(DP), INTENT(OUT) :: StopTime !<On return, the stop time for the time control loop.
    REAL(DP), INTENT(OUT) :: TimeIncrement !<On return, the time increment for the time control loop.
    REAL(DP), INTENT(OUT) :: CurrentTime !<On return, the current time for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTimesGetPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_TIMES_GET(ControlLoop,StartTime,StopTime,TimeIncrement,CurrentTime,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTimesGetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesGetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesGetPtr

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

    CALL ENTERS("CMISSControlLoopTimesSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,Err,ERROR,*999)

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
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("CMISSControlLoopTimesSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,StartTime,StopTime,TimeIncrement,Err,ERROR,*999)

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
  
  !>Sets/changes the time parameters for a time control loop identified by a pointer.
  SUBROUTINE CMISSControlLoopTimesSetPtr(ControlLoop,StartTime,StopTime,TimeIncrement,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to set the times for.
    REAL(DP), INTENT(IN) :: StartTime !<The start time for the time control loop to set.
    REAL(DP), INTENT(IN) :: StopTime !<The stop time for the time control loop to set.
    REAL(DP), INTENT(IN) :: TimeIncrement !<The time increment for the time control loop to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTimesSetPtr",Err,ERROR,*999)
    
    CALL CONTROL_LOOP_TIMES_SET(ControlLoop,StartTime,StopTime,TimeIncrement,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTimesSetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopTimesSetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopTimesSetPtr")
    CALL CMISS_HANDLE_ERRORS(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTimesSetPtr

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

    CALL ENTERS("CMISSControlLoopTypeSetNumber0",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifier,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,LoopType,Err,ERROR,*999)

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

    CALL ENTERS("CMISSControlLoopTypeSetNumber1",Err,ERROR,*999)
    
    NULLIFY(CONTROL_LOOP)
    NULLIFY(PROBLEM)
    CALL PROBLEM_USER_NUMBER_FIND(ProblemUserNumber,PROBLEM,Err,ERROR,*999)
    CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,ControlLoopIdentifiers,CONTROL_LOOP,Err,ERROR,*999)
    CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,LoopType,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTypeSetNumber1")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeSetNumber1",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeSetNumber1")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeSetNumber0

  !  
  !================================================================================================================================
  !  
  
  !>Sets/changes the loop type for a control loop identified by a pointer. \todo is this really public???
  SUBROUTINE CMISSControlLoopTypeSetPtr(ControlLoop,LoopType,Err)
  
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), INTENT(IN), POINTER :: ControlLoop !<A pointer to the control loop to set the loop type for.
    INTEGER(INTG), INTENT(IN) :: LoopType !<The type of control loop to set. \see OPENCMISS_ProblemControlLoopTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSControlLoopTypeSetPtr",Err,ERROR,*999)

    CALL CONTROL_LOOP_TYPE_SET(ControlLoop,LoopType,Err,ERROR,*999)

    CALL EXITS("CMISSControlLoopTypeSetPtr")
    RETURN
999 CALL ERRORS("CMISSControlLoopTypeSetPtr",Err,ERROR)
    CALL EXITS("CMISSControlLoopTypeSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSControlLoopTypeSetPtr

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

    CALL ENTERS("CMISSCoordinateSystemCreateFinishNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,Err,ERROR,*999)

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
  
  !>Finishes the creation of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemCreateFinishPtr(CoordinateSystem,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN), POINTER :: CoordinateSystem !<A pointer to the coordinate system to finish creating.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSCoordinateSystemCreateFinishPtr",Err,ERROR,*999)
 
    CALL COORDINATE_SYSTEM_CREATE_FINISH(CoordinateSystem,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemCreateFinishPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemCreateFinishPtr",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemCreateFinishPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemCreateFinishPtr

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
  
  !>Starts the creation of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemCreateStartPtr(CoordinateSystemUserNumber,CoordinateSystem,Err)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemUserNumber !<The user number of the coordinate system to start creating.
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(OUT), POINTER :: CoordinateSystem !<On return, a pointer to the coordinate system that has been created.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSCoordinateSystemCreateStartPtr",Err,ERROR,*999)
 
    CALL COORDINATE_SYSTEM_CREATE_START(CoordinateSystemUserNumber,CoordinateSystem,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemCreateStartPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemCreateStartPtr",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemCreateStartPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemCreateStartPtr

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

    CALL ENTERS("CMISSCoordinateSystemDestroyNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_DESTROY(COORDINATE_SYSTEM,Err,ERROR,*999)

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
  
  !>Destroys a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemDestroyPtr(CoordinateSystem,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(INOUT), POINTER :: CoordinateSystem !<A pointer to the coordinate system to destroy.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables
 
    CALL ENTERS("CMISSCoordinateSysteDestroyPtr",Err,ERROR,*999)
 
    CALL COORDINATE_SYSTEM_DESTROY(CoordinateSystem,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemDestroyPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDestroyPtr",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDestroyPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDestroyPtr

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

    CALL ENTERS("CMISSCoordinateSystemDimensionGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,CoordinateSystemDimension,Err,ERROR,*999)

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
 
  !>Returns the dimension of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemDimensionGetPtr(CoordinateSystem,CoordinateSystemDimension,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to get the dimension for.
    INTEGER(INTG), INTENT(OUT) :: CoordinateSystemDimension !<On return, the dimension of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemDimensionGetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_DIMENSION_GET(CoordinateSystem,CoordinateSystemDimension,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemDimensionGetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDimensionGetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemDimensionGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDimensionGetPtr

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

    CALL ENTERS("CMISSCoordinateSystemDimensionSetNumber",Err,ERROR,*999)
    
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,CoordinateSystemDimension,Err,ERROR,*999)

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
 
  !>Sets/changes the dimension of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemDimensionSetPtr(CoordinateSystem,CoordinateSystemDimension,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to set the dimension for.
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemDimension !<The dimension of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemDimensionSetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_DIMENSION_SET(CoordinateSystem,CoordinateSystemDimension,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemDimensionSetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemDimensionSetPtr",Err,ERROR)
    CALL EXITS("CMISSCoordinateSystemDimensionSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemDimensionSetPtr

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

    CALL ENTERS("CMISSCoordinateSystemFocusGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_FOCUS_GET(COORDINATE_SYSTEM,Focus,Err,ERROR,*999)

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
 
  !>Returns the focus of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemFocusGetPtr(CoordinateSystem,Focus,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to get the focus for.
    REAL(DP), INTENT(OUT) :: Focus !<On return, the focus of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemFocusGetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_FOCUS_GET(CoordinateSystem,Focus,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemFocusGetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemFocusGetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemFocusGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemFocusGetPtr

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

    CALL ENTERS("CMISSCoordinateSystemFocusSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_FOCUS_SET(COORDINATE_SYSTEM,Focus,Err,ERROR,*999)

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
 
  !>Sets/changes the focus of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemFocusSetPtr(CoordinateSystem,Focus,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to set the focus for.
    REAL(DP), INTENT(IN) :: Focus !<The focus of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemFocusSetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_FOCUS_SET(CoordinateSystem,Focus,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemFocusSetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemFocusSetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemFocusSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemFocusSetPtr

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

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_GET(COORDINATE_SYSTEM,RadialInterpolationType,Err,ERROR,*999)

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
 
  !>Returns the radial interpolation type of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemRadialInterpolationGetPtr(CoordinateSystem,RadialInterpolationType,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to get the radial interpolation type for.
    INTEGER(INTG), INTENT(OUT) :: RadialInterpolationType !<On return, the radial interpolation type of the coordinate system. \see OPENCMISS_CoordinateRadialInterpolations
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationGetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_GET(CoordinateSystem,RadialInterpolationType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemRadialInterpolationGetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemRadialInterpolationGetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemRadialInterpolationGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemRadialInterpolationGetPtr

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

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_SET(COORDINATE_SYSTEM,RadialInterpolationType,Err,ERROR,*999)

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
 
  !>Sets/changes the radial interpolation type of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemRadialInterpolationSetPtr(CoordinateSystem,RadialInterpolationType,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to set the radial interpolation type for.
    INTEGER(INTG), INTENT(IN) :: RadialInterpolationType !<The radial interpolation type of the coordinate system to set. \see OPENCMISS_CoordinateRadialInterpolations
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemRadialInterpolationSetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_SET(CoordinateSystem,RadialInterpolationType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemRadialInterpolationSetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemRadialInterpolationSetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemRadialInterpolationSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemRadialInterpolationSetPtr

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

    CALL ENTERS("CMISSCoordinateSystemTypeGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_TYPE_GET(COORDINATE_SYSTEM,CoordinateSystemType,Err,ERROR,*999)

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
 
  !>Returns the type of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemTypeGetPtr(CoordinateSystem,CoordinateSystemType,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to get the type for.
    INTEGER(INTG), INTENT(OUT) :: CoordinateSystemType !<On return, the type of the coordinate system. \see OPENCMISS_CoordinateSystemTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemTypeGetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_TYPE_GET(CoordinateSystem,CoordinateSystemType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemTypeGetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeGetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemTypeGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeGetPtr

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

    CALL ENTERS("CMISSCoordinateSystemTypeSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_TYPE_SET(COORDINATE_SYSTEM,CoordinateSystemType,Err,ERROR,*999)

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
 
  !>Sets/changes the type of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemTypeSetPtr(CoordinateSystem,CoordinateSystemType,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to set the type for.
    INTEGER(INTG), INTENT(IN) :: CoordinateSystemType !<The type of the coordinate system to set. \see OPENCMISS_CoordinateSystemTypes
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemTypeSetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_TYPE_SET(CoordinateSystem,CoordinateSystemType,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemTypeSetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemTypeSetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemTypeSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemTypeSetPtr
  
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

    CALL ENTERS("CMISSCoordinateSystemOriginGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_ORIGIN_GET(COORDINATE_SYSTEM,Origin,Err,ERROR,*999)

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
 
  !>Returns the origin of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemOriginGetPtr(CoordinateSystem,Origin,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to get the origin for.
    REAL(DP), INTENT(OUT) :: Origin(:) !<On return, the origin of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOriginGetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIGIN_GET(CoordinateSystem,Origin,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOriginGetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOriginGetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOriginGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOriginGetPtr

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

    CALL ENTERS("CMISSCoordinateSystemOriginSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_ORIGIN_SET(COORDINATE_SYSTEM,Origin,Err,ERROR,*999)

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
 
  !>Sets/changes the origin of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemOriginSetPtr(CoordinateSystem,Origin,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to set the origin for.
    REAL(DP), INTENT(IN) :: Origin(:) !<The origin of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOriginSetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIGIN_SET(CoordinateSystem,Origin,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOriginSetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOriginSetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOriginSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOriginSetPtr

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

    CALL ENTERS("CMISSCoordinateSystemOrientationGetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_ORIENTATION_GET(COORDINATE_SYSTEM,ORIENTATION,Err,ERROR,*999)

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
 
  !>Returns the orientation of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemOrientationGetPtr(CoordinateSystem,Orientation,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to get the orientation for.
    REAL(DP), INTENT(OUT) :: Orientation(:,:) !<On return, the orientation of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOrientationGetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIENTATION_GET(CoordinateSystem,Orientation,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOrientationGetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOrientationGetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOrientationGetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOrientationGetPtr

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

    CALL ENTERS("CMISSCoordinateSystemOrientationSetNumber",Err,ERROR,*999)
 
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(CoordinateSystemUserNumber,COORDINATE_SYSTEM,Err,ERROR,*999)
    CALL COORDINATE_SYSTEM_ORIENTATION_SET(COORDINATE_SYSTEM,ORIENTATION,Err,ERROR,*999)

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
 
  !>Sets/changes the orientation of a coordinate system identified by a pointer.
  SUBROUTINE CMISSCoordinateSystemOrientationSetPtr(CoordinateSystem,Orientation,Err)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM), INTENT(IN), POINTER  :: CoordinateSystem !<A pointer to the coordinate system to set the orientation for.
    REAL(DP), INTENT(IN) :: Orientation(:,:) !<The orientation of the coordinate system to set.
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
    !Local variables

    CALL ENTERS("CMISSCoordinateSystemOrientationSetPtr",Err,ERROR,*999)
    
    CALL COORDINATE_SYSTEM_ORIENTATION_SET(CoordinateSystem,Orientation,Err,ERROR,*999)

    CALL EXITS("CMISSCoordinateSystemOrientationSetPtr")
    RETURN
999 CALL ERRORS("CMISSCoordinateSystemOrientationSetPtr",Err,ERROR)
    CALL EXIT("CMISSCoordinateSystemOrientationSetPtr")
    CALL CMISS_HANDLE_ERROR(Err,ERROR)
    RETURN
    
  END SUBROUTINE CMISSCoordinateSystemOrientationSetPtr

!!==================================================================================================================================
!!
!! FIELD_ROUTINES
!!
!!==================================================================================================================================

  
 
END MODULE OPENCMISS
