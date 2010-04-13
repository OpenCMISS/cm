/* 
 * \file
 * $Id: opencmiss.h 582 2009-07-08 03:51:19Z catalept $
 * \author Chris Bradley
 * \brief The OpenCMISS library C header file.
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is OpenCMISS
 *
 * The Initial Developer of the Original Code is University of Auckland,
 * Auckland, New Zealand and University of Oxford, Oxford, United
 * Kingdom. Portions created by the University of Auckland and University
 * of Oxford are Copyright (C) 2007 by the University of Auckland and
 * the University of Oxford. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */

#ifndef OPENCMISS_H
#define OPENCMISS_H
#endif

/* 
 * Defines
 */
const int CMISSNoError=0;

const int CMISSTrue=1;
const int CMISSFalse=0;

/*
 *==================================================================================================================================
 *
 * ANALYTIC_ANALYSIS_ROUTINES
 *
 *==================================================================================================================================
 */

/*
 *==================================================================================================================================
 *
 * BASE_ROUTINES
 *
 *==================================================================================================================================
 */

/* Parameters */

/* \addtogroup OPENCMISS_DiagnosticAndTimingConstants OPENCMISS::DiagnosticAndTiming::Constants
 * \brief Diagnostic and Timing constants.
 *@{
 * \addtogroup OPENCMISS_DiagnosticTypes OPENCMISS::DiagnosticAndTiming::DiagnosticTypes
 * \brief Diganostic constants.
 * \see OPENCMISS::DiagnosticTypes,OPENCMISS
 *@{
 */
#define CMISSAllDiagType 1;
#define CMISSInDiagType 2;
#define CMISSFromDiagType 3;
/* @}
 * \addtogroup OPENCMISS_TimingTypes OPENCMISS::DiagnosticAndTiming::TimingTypes
 * \brief Timing constants.
 * \see OPENCMISS::TimingTypes,OPENCMISS
 *@{
 */
static int CMISSAllTimingType = 1;
static int CMISSInTimingType = 2;
static int CMISSFromTimingType = 3;
/*@}
 @}*/

/*
 *==================================================================================================================================
 *
 * BASIS_ROUTINES
 *
 *==================================================================================================================================
 */

/* Module parameters */

/* > \addtogroup OPENCMISS_BasisConstants OPENCMISS::Basis::Constants
 *> \brief Basis function constants.
 *>@{
 *> \addtogroup OPENCMISS_BasisTypes OPENCMISS::Basis::BasisTypes
 *> \brief Basis definition type parameters.
 *> \see OPENCMISS::BasisConstants,OPENCMISS
 *>@{
 */
static int CMISSBasisLagrangeHermiteTPType = 1;
static int CMISSBasisSimplexType = 2;
static int CMISSBasisSerendipityType = 3;
static int CMISSBasisAuxilliaryType = 4;
static int CMISSBSplineTPType = 5;
static int CMISSBasisFourierLagrangeHermiteTPType = 6;
static int CMISSBasisExtendedLagrangeTPType = 7;
/*
 *>@}
 *> \addtogroup OPENCMISS_BasisInterpolationSpecifications OPENCMISS::Basis::InterpolationSpecifications
 *> \brief Interpolation specification parameters
 *> \see OPENCMISS::BasisConstants,OPENCMISS
 *>@{
 */
static int CMISSBasisLinearLagrangeInterpolation = 1;
static int CMISSBasisQuadraticLagrangeInterpolation = 2;
static int CMISSBasisCubicLagrangeInterpolation = 3;
static int CMISSBasisCubicHermiteInterpolation = 4;
static int CMISSBasisQuadratic1HermiteInterpolation = 5;
static int CMISSBasisQuadratic2HermiteInterpolation = 6;
static int CMISSBasisLinearSimplexInterpolation = 7;
static int CMISSBasisQuadraticSimplexInterpolation = 8;
static int CMISSBasisCubicSimplexInterpolation = 9;
/*
 * >@}
 * > \addtogroup OPENCMISS_BasisQuadratureTypes OPENCMISS::Basis::QuadratureTypes
 * > \brief Basis quadrature type parameters.
 * > \see OPENCMISS::BasisConstants,OPENCMISS
 * >@{
 */
static int CMISSBasisGaussLegendreQuadrature = 1;
static int CMISSBasisGaussLaguerreQuadrature = 2;
static int CMISSBasisGaussHermiteQuadrature = 3;
static int CMISSBasisAdaptiveGaussLegendreQuadrature = 4;
static int CMISSBasisGaussSimplexQuadratire = 5;
/*
 * >@}
 * > \addtogroup OPENCMISS_BasisXiCollapse OPENCMISS::Basis::XiCollapse
 * > \brief Basis Xi collapse parameters.
 * > \see OPENCMISS::Basis,OPENCMISS
 * >@{
 */
static int CMISSBasisXiCollapsed = 1;
static int CMISSBasisCollapsedAtXi0 = 2;
static int CMISSBasisCollapsedAtXi1 = 3;
static int CMISSBasisNotCollapsed = 4;
/*
 * >@}
 * >@}
 */

/*
 *==================================================================================================================================
 *
 * BOUNDARY_CONDITIONS_ROUTINES
 *
 *==================================================================================================================================
 */

/* Module parameters */


/*
 * > \addtogroup OPENCMISS_BoundaryConditionsConstants OPENCMISS::BoundaryConditions::Constants
 * > \brief Boundary conditions constants.
 * >@{
 * > \addtogroup OPENCMISS_BoundaryConditionsTypes OPENCMISS::BoundaryConditions::Types
 * > \brief Boundary conditions type parameters.
 * > \see OPENCMISS::BoundaryConditions,OPENCMISS
 * >@{
 */
static int CMISSBoundaryConditionNotFixed = 0;
static int CMISSBoundaryConditionFixed = 1;
static int CMISSBoundaryConditionMixed = 7;
/* Temporary boundary flags (to be removed when general boundary object becomes available!) */
static int CMISSBoundaryConditionFixedWall = 4;
static int CMISSBoundaryConditionInletWall = 2;
static int CMISSBoundaryConditionMovedWall = 5;
/*
 *>@}
 *>@}
 */

/*
 *==================================================================================================================================
 *
 * CMISS
 *
 *==================================================================================================================================
 */

static int CMISSReturnErrorCode = 0;
static int CMISSOutputError = 1;
static int CMISSTrapError = 2;

/*
 *==================================================================================================================================
 *
 * COMP_ENVIRONMENT
 *
 *==================================================================================================================================
 */

/*
 *==================================================================================================================================
 *
 * CONSTANTS
 *
 *==================================================================================================================================
 */

/* Module parameters */

/* > \addtogroup OPENCMISS_Constants OPENCMISS::Constants
 * > \brief Control loops constants.
 * >@{
 * > \addtogroup OPENCMISS_GlobalDerivativeConstants OPENCMISS::Constants::GlobalDerivativeConstants
 * > \brief Global derivative constant identifiers
 * > \see OPENCMISS_CONSTANTS,OPENCMISS
 * >@{
 */
static int CMISSNoGlobalDerivative = 1;
static int CMISSGlobalDerivativesS1 = 2;
static int CMISSGlobalDerivativesS2 = 3;
static int CMISSGlobalDerivativesS1S2 = 4;
static int CMISSGlobalDerivativesS3 = 5;
static int CMISSGlobalDerivativesS1S3 = 6;
static int CMISSGlobalDerivativesS2S3 = 7;
static int CMISSGlobalDerivativesS1S2S3 = 8;
/*
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 * CONTROL_LOOP_ROUTINES
 *
 *==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_ControlLoopConstants OPENCMISS::ControlLoop::Constants
 * > \brief Control loops constants.
 * >@{
 * > \addtogroup OPENCMISS_ControlLoopIdentifiers OPENCMISS::ControlLoop::Identifiers
 * > \brief The control loop identification parameters.
 * > \see OPENCMISS::ControlLoop,OPENCMISS
 * >@{
 */
static int CMISSControlLoopNode = 0;

/*
 * >@}
 * >@}
 */

/* ==================================================================================================================================
 *
 * COORDINATE_ROUTINES
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_CoordinateConstants OPENCMISS::Coordinate::Constants
 * > \brief Coordinate constants.
 * >@{
 * > \addtogroup OPENCMISS_CoordinateSystemTypes OPENCMISS::Coordinate::SystemTypes
 * > \brief Coordinate system type parameters.
 * > \see OPENCMISS::Coordinate,OPENCMISS
 * >@{
 */
static int CMISSCoordinateRectangularCartesianType = 1;
static int CMISSCoordinateCylindricalPolarType = 2;
static int CMISSCoordinateSphericalPolarType = 3;
static int CMISSCoordinateProlateSpheroidalType = 4;
static int CMISSCoordinateOblateSpheroidalType = 5;
/*
 * >@}
 * > \addtogroup OPENCMISS_CoordinateRadialInterpolations OPENCMISS::Coordinate::RadialInterpolations
 * > \brief The type of radial interpolation for polar coordinate systems
 * > \see OPENCMISS::Coordinate,OPENCMISS
 * >@{
 */
static int CMISSCoordinateNoRadialInterpolationType = 0;
static int CMISSCoordinateRadialInterpolationType = 1;
static int CMISSCoordinateRadialSquaredInterpolationType = 2;
static int CMISSCoordinateRadialCubedInterpolationType = 3;

/*
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 * EQUATIONS_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_EquationsConstants OPENCMISS::Equations::Constants
 * > \brief Equations  constants.
 * >@{
 * > \addtogroup OPENCMISS_EquationsOutputTypes OPENCMISS::Equations::OutputTypes
 * > \brief Equations output types
 * > \see OPENCMISS::Equations,OPENCMISS
 * >@{
 */
static int CMISSEquationsNoOutput = 0;
static int CMISSEquationsTimingOutput = 1;
static int CMISSEquationsMatrixOutput = 2;
static int CMISSEquationsElementMatrixOutput = 3;

/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsSparsityTypes OPENCMISS::Equations::SparsityTypes
 * > \brief Equations sparsity types
 * > \see OPENCMISS::Equations,OPENCMISS
 * >@{
 */
static int CMISSEquationsSparseMatrices = 1;
static int CMISSEquationsFullMatrices = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsLumpingTypes OPENCMISS::Equations::LumpingTypes
 * > \brief Equations lumping types
 * > \see OPENCMISS::Equations,OPENCMISS
 * >@{
 */
static int CMISSEquationsUnlumpedMatrices = 1;
static int CMISSEquationsLumpedMatrices = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsLinearityTypes OPENCMISS::Equations::LinearityTypes
 * > \brief The equations linearity types
 * > \see OPENCMISS::Equations,OPENCMISS
 * >@{
 */
static int CMISSEquationsLinear = 1;
static int CMISSEquationsNonlinear = 2;
static int CMISSEquationsNonlinearBCs = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsTimeDepedenceTypes OPENCMISS::Equations::TimeDepedenceTypes
 * > \brief The equations time dependence types
 * > \see OPENCMISS::Equations,OPENCMISS
 * >@{
 */
static int CMISSEquationsStatic = 1;
static int CMISSEquationsQuasistatic = 2;
static int CMISSEquationsFirstOrderDynamic = 3;
static int CMISEquationsSecondOrderDynamic = 4;
static int CMISSEquationsTimeStepping = 5;
/*
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 *
 * EQUATIONS_SET_CONSTANTS
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_EquationsSetConstants OPENCMISS::EquationsSet::Constants
 * > \brief Equations set constants.
 * >@{
 * > \addtogroup OPENCMISS_EquationsSetClasses OPENCMISS::EquationsSet::Classes
 * > \brief Equations set classes.
 * > \see OPENCMISS::EquationsSet,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetNoClass = 0;
static int CMISSEquationsSetElasticityClass = 1;
static int CMISSEquationsSetFluidMechanicsClass = 2;
static int CMISSEquationsSetElectroMechanicsClass = 3;
static int CMISSEquationsSetClassicalFieldClass = 4;
static int CMISSEquationsSetBioelectricsClass = 5;
static int CMISSEquationsSetModalClass = 6;
static int CMISSEquationsSetFittingClass = 7;
static int CMISSEquationsSetOptimisationClass = 8;
/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsSetTypes OPENCMISS::EquationsSet::Types
 * > \brief Equations set Types.
 * > \see OPENCMISS::EquationsSet,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetNoType = 0;
static int CMISSEquationsSetLinearElasticityType = 1;
static int CMISSEquationsSetFiniteElasticityType = 2;
static int CMISSEquationsSetStokesEquationType = 1;
static int CMISSEquationsSetNavierStokesEquationType = 2;
static int CMISSEquationsSetDarcyEquationType = 3;
static int CMISSEquationsSetElectrostaticType = 1;
static int CMISSEquationsSetMagnetostaticType = 2;
static int CMISSEquationsSetMaxwellsEquationType = 3;
static int CMISSEquationsSetLaplaceEquationType = 1;
static int CMISSEquationsSetPoissonEquationType = 2;
static int CMISSEquationsSetHelmholtzEquationType = 3;
static int CMISSEquationsSetWaveEquationType = 4;
static int CMISSEquationsSetDiffusionEquationType = 5;
static int CMISSEquationsSetAdvectionDiffusionEquationType = 6;
static int CMISSEquationsSetReactionDiffusionEquationType = 7;
static int CMISSEquationsSetBiharmonicEquationType = 8;
static int CMISSEquationsSetMonodomainEquationType = 1;
static int CMISSEquationsSetBidomainEquationType = 2;
static int CMISSEquationsSetLinearElasticModalType = 1;
static int CMISSEquationsSetGalerkinProjectionEquationType = 9;
/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsSetSubtypes OPENCMISS::EquationsSet::Subtypes
 * > \brief Equations set subtypes.
 * > \see OPENCMISS::EquationsSet,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetNoSubtype = 0;
static int CMISSEquationsSetThreeDimensionalSubtype = 1;
static int CMISSEquationsSetPlaneStressSubtype = 2;
static int CMISSEquationsSetPlaneStrainSubtype = 3;
static int CMISSEquationsSetOneDimensionalSubtype = 4;
static int CMISSEquationsSetPlateSubtype = 5;
static int CMISSEquationsSetShellSubtype = 6;
static int CMISSEquationsSetStaticStokesSubtype = 1;
static int CMISSEquationsSetLaplaceStokesSubtype = 1;
static int CMISSEquationsSetTransientStokesSubtype = 3;
static int CMISSEquationsSetALEStokesSubtype = 5;
static int CMISSEquationsSetOptimisedStokesSubtype = 4;
static int CMISSEquationsSetStaticNavierStokesSubtype = 1;
static int CMISSEquationsSetLaplaceNavierStokesSubtype = 2;
static int CMISSEquationsSetTransientNavierStokesSubtype = 3;
static int CMISSEquationsSetALENavierStokesSubtype = 5;
static int CMISSEquationsSetOptimisedNavierStokesSubtype = 4;
static int CMISSEquationsSetStandardDarcySubtype = 1;
static int CMISSEquationsSetQuasistaticDarcySubtype = 2;
static int CMISSEquationsSetALEDarcySubtype = 3;
static int CMISSEquationsSetStandardLaplaceSubtype = 1;
static int CMISSEquationsSetGeneralisedLaplaceSubtype = 2;
static int CMISSEquationsSetMovingMeshLaplaceSubtype = 3;
static int CMISSEquationsSetConstantSourcePoissonSubtype = 1;
static int CMISSEquationsSetLinearSourcePoissonSubtype = 2;
static int CMISSEquationsSetQuadraticSourcePoissonSubtype = 3;
static int CMISSEquationsSetExponentialSourcePoissonSubtype = 4;
static int CMISSEquationsSetStokesPoissonSubtype = 5;
static int CMISSEquationsSetNavierStokesPoissonSubtype = 6;
static int CMISSEquationsSetNoSourceHelmholtzSubtype = 1;
static int CMISSEquationsSetNoSourceDiffusionSubtype = 1;
static int CMISSEquationsSetConstantSourceDiffusionSubtype = 2;
static int CMISSEquationsSetLinearSourceDiffusionSubtype = 3;
static int CMISSEquationsSetQuadraticSourceDiffusionSubtype = 4;
static int CMISSEquationsSetExponentialSourceDiffusionSubtype = 5;
static int CMISSEquationsSetNoSourceAdvectionDiffusionSubtype = 1;
static int CMISSEquationsSetConstantSourceAdvectionDiffusionSubtype = 2;
static int CMISSEquationsSetLinearSourceAdvectionDiffusionSubtype = 3;
static int CMISSEquationsSetQuadraticSourceAdvectionDiffusionSubtype = 4;
static int CMISSEquationsSetExponentialSourceAdvectionDiffusionSubtype = 5;
static int CMISSEquationsSetFirstBidomainSubtype = 1;
static int CMISSEquationsSetSecondBidomainSubtype = 2;
static int CMISSEquationsSetStandardGalerkinProjectionSubtype = 1;
static int CMISSEquationsSetGeneralisedGalerkinProjectionSubtype = 2;
static int CMISSEquationsSetMatPropertiesGalerkinProjectionSubtype = 3;

/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsSetSolutionMethods OPENCMISS::EquationsSet::SolutionMethods
 * > \brief The solution method parameters
 * > \see OPENCMISS::EquationsSet,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetFEMSolutionMethod = 1;
static int CMISSEquationsSetBEMSolutionMethod = 2;
static int CMISSEquationsSetFDSolutionMethod = 3;
static int CMISSEquationsSetFVSolutionMethod = 4;
static int CMISSEquationsSetGFEMSolutionMethod = 5;
static int CMISSEquationsSetGFDSolutionMethod = 6;
static int CMISSEquationsSetGFVSolutionMethod = 7;
/*
 * >@}
 * > \addtogroup OPENCMISS_EquationsSetAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes
 * > \brief The analytic function types.
 * > \see OPENCMISS::EquationsSet,OPENCMISS
 * >@{
 * > \addtogroup OPENCMISS_EquationsSetLaplaceAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::Laplace
 * > \brief The analytic function types for a Laplace equation
 * > \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetLaplaceEquationTwoDim1 = 1;
static int CMISSEquationsSetLaplaceEquationTwoDim2 = 2;
static int CMISSEquationsSetLaplaceEquationThreeDim1 = 3;
static int CMISSEquationsSetLaplaceEquationThreeDim2 = 4;
/*
 * >@}
 * > \addtogroup OPENCMISS_PoissonAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::Poisson
 * > \brief The analytic function types for a Poisson equation.
 * > \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetPoissonTwoDim1 = 1;
static int CMISSEquationsSetPoissonTwoDim2 = 2;
static int CMISSEquationsSetPoissonTwoDim3 = 3;
static int CMISSEquationsSetPoissonThreeDim1 = 4;
static int CMISSEquationsSetPoissonThreeDim2 = 5;
static int CMISSEquationsSetPoissonThreeDim3 = 6;
/*
 * >@}
 * > \addtogroup OPENCMISS_PoissonAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::Poisson
 * > \brief The analytic function types for a Poisson equation.
 * > \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetDiffusionTwoDim1 = 1;
static int CMISSEquationsSetLinearSourceDiffusionThreeDim1 = 2;
/*
 * > \addtogroup OPENCMISS_StokesAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::Stokes
 * > \brief The analytic function types for a Stokes equation.
 * > \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetStokesTwoDim1 = 1;
static int CMISSEquationsSetStokesTwoDim2 = 2;
static int CMISSEquationsSetStokesTwoDim3 = 3;
static int CMISSEquationsSetStokesThreeDim1 = 4;
static int CMISSEquationsSetStokesThreeDim2 = 5;
static int CMISSEquationsSetStokesThreeDim3 = 6;
/*
 * >@}
 * > \addtogroup OPENCMISS_NavierStokesAnalyticFunctionTypes OPENCMISS::EquationsSet::AnalyticFunctionTypes::NavierStokes
 * > \brief The analytic function types for a Navier-Stokes equation.
 * > \see OPENCMISS::EquationsSet::AnalyticFunctionTypes,OPENCMISS
 * >@{
 */
static int CMISSEquationsSetNavierStokesTwoDim1 = 1;
static int CMISSEquationsSetNavierStokesTwoDim2 = 2;
static int CMISSEquationsSetNavierStokesTwoDim3 = 3;
static int CMISSEquationsSetNavierStokesThreeDim1 = 4;
static int CMISSEquationsSetNavierStokesThreeDim2 = 5;
static int CMISSEquationsSetNavierStokesThreeDim3 = 6;
static int CMISSEquationsSetNavierStokesTwoDim4 = 7;
/* >@}
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 * EQUATIONS_SET_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */


/*
 * ==================================================================================================================================
 *
 *  FIELD_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_FieldConstants OPENCMISS::Field::Constants
 * > \brief Field constants.
 * >@{
 * > \addtogroup OPENCMISS_FieldDependentTypes OPENCMISS::Field::DependentTypes
 * > \brief Depedent field parameter types.
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldIndependentType = 1;
static int CMISSFieldDependentType = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldDimensionTypes OPENCMISS::Field::DimensionTypes
 * > \brief Field dimension parameter types.
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldScalarDimensionType = 1;
static int CMISSFieldVectorDimensionType = 2;
static int CMISSFieldTensorDimensionType = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldTypes OPENCMISS::Field::Types
 * > \brief Field type parameters.
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldGeometricType = 1;
static int CMISSFieldFibreType = 2;
static int CMISSFieldGeneralType = 3;
static int CMISSFieldMaterialType = 4;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldInterpolationTypes OPENCMISS::Field::InterpolationTypes
 * > \brief Field interpolation parameters.
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldConstantInterpolation = 1;
static int CMISSFieldElementBasedInterpolation = 2;
static int CMISSFieldNodeBasedInterpolation = 3;
static int CMISSFieldGridPointBasedInterpolation = 4;
static int CMISSFieldGaussPointBasedInterpolation = 5;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldVariableTypes OPENCMISS::Field::VariableTypes
 * > \brief Field variable type parameters.
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldUVariableType = 1;
static int CMISSFieldDelUDelNVariableType = 2;
static int CMISSFieldDelUDelTVariableType = 3;
static int CMISSFieldDel2UDelT2VariableType = 4;
static int CMISSFieldVVariableType = 5;
static int CMISSFieldDelVDelNVariableType = 6;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldDataTypes OPENCMISS::Field::DataTypes
 * > \brief Field data types
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldIntgType = 1;
static int CMISSFieldSPType = 2;
static int CMISSFieldDPType = 3;
static int CMISSFieldLType = 4;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldDOFOrderTypes OPENCMISS::Field::DOFOrderTypes
 * > \brief Field DOF order types
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldSeparatedComponentDOFOrder = 1;
static int CMISSFieldContiguousComponentDOFOrder = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldParameterSetTypes OPENCMISS::Field::ParameterSetTypes
 * > \brief Field parameter set type parameters
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISSFieldValuesSetType = 1;
static int CMISSInitialValuesSetType = 3;
static int CMISSFieldIncrementalValuesSetType = 4;
static int CMISSFieldAnalyticValuesSetType = 5;
static int CMISSPreviousValuesSetType = 6;
static int CMISSMeanPredictedDisplacementSetType = 7;
static int CMISSFieldVelocityValuesSetType = 8;
static int CMISSFieldInitialVelocitySetType = 9;
static int CMISSFieldPreviousVelocitySetType = 9;
static int CMISSFieldMeanPredictedVelocitySetType = 10;
static int CMISSFieldAccelerationValuesSetType = 11;
static int CMISSInitialAccelerationSetType = 12;
static int CMISSFieldPreviousAccelerationSetType = 12;
static int CMISSMeanPredictedAccelerationSetType = 13;
/*
 * >@}
 * > \addtogroup OPENCMISS_FieldScalingTypes OPENCMISS::Field::ScalingTypes
 * > \brief Field scaling type parameters
 * > \see OPENCMISS::Field,OPENCMISS
 * >@{
 */
static int CMISFieldNoScaling = 0;
static int CMISSFieldUnitScaling = 1;
static int CMISSFieldArcLengthScaling = 2;
static int CMISSFieldArithmeticMeanScaling = 3;
static int CMISSFieldHarmonicMeanScaling = 4;
/*
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 *  FIELD_IO_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * ==================================================================================================================================
 *
 * GENERATED_MESH_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_GeneratedMeshConstants OPENCMISS::GeneratedMesh::Constants
 * > \brief Generated mesh constants.
 * >@{
 * > \addtogroup OPENCMISS_GeneratedMeshTypes OPENCMISS::GeneratedMesh::Types
 * > \brief Generated mesh types.
 * > \see OPENCMISS::GeneratedMesh,OPENCMISS
 * >@{
 */
static int CMISSGeneratedMeshRegularMeshType = 1;
static int CMISSGeneratedMeshPolarMeshType = 2;
static int CMISSGeneratedMeshFractalTreeMeshType = 3;
/*
 * >@}
 * >@}
 */

/* ==================================================================================================================================
 *
 * KINDS
 *
 * ==================================================================================================================================
 */



/*
 * ==================================================================================================================================
 *
 *  MESH_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_MeshConstants OPENCMISS::Mesh::Constants
 * > \brief Mesh constants.
 * >@{
 * > \addtogroup OPENCMISS_DecompositionTypes OPENCMISS::Mesh::DecompositionTypes
 * > \brief The Decomposition types parameters
 * > \see OPENCMISS::Mesh,OPENCMISS
 * >@{
 */
static int CMISSDecompositionAllType = 1;
static int CMISSDecompositionCalculatedType = 2;
static int CMISSDecompositionUserDefinedType = 3;
/*
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 *  NODE_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */


/*
 * ==================================================================================================================================
 *
 *  PROBLEM_CONSTANTS_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_ProblemConstants OPENCMISS::Problem::Constants
 * > \brief Problem constants.
 * >@{
 * > \addtogroup OPENCMISS_ProblemClasses OPENCMISS::Problem::Classes
 * > \brief Problem classes.
 * > \see OPENCMISS::Problem,OPENCMISS
 * >@{
 */
static int CMISSProblemNoClass = 0;
static int CMISSProbelmElasticityClass = 1;
static int CMISSProblemFluidMechanicsClass = 2;
static int CMISSProblemElectromagneticsClass = 3;
static int CMISSProblemClassicalFieldClass = 4;
static int CMISSProblemBioelectricsClass = 5;
static int CMISSProblemModalClass = 6;
static int CMISSProblemFittingClass = 7;
static int CMISSProblemOptimisationClass = 8;
static int CMISSProblemMultiPhysicsClass=9;
/*
 * >@}
 * > \addtogroup OPENCMISS_ProblemTypes OPENCMISS::Problem::Types
 * > \brief Problem Types.
 * > \see OPENCMISS::Problem,OPENCMISS
 * >@{
 */
static int CMISSProblemNoType = 0;
static int CMISSProblemLinearElasticityType = 1;
static int CMISSProblemFiniteElasticityType = 2;
static int CMISSProblemStokesEquationType = 1;
static int CMISSProblemNavierStokesEquationType = 2;
static int CMISSProblemDarcyEquationType = 3;
static int CMISSProblemElectrostaticType = 1;
static int CMISSProblemMagnetostaticType = 2;
static int CMISSProblemMaxwellsEquationsType = 3;
static int CMISSProblemLaplaceEquationType = 1;
static int CMISSProblemPoissonEquationType = 2;
static int CMISSProblemHelmholtzEquationType = 3;
static int CMISSProblemWaveEquationType = 4;
static int CMISSProblemDiffusionEquationType = 5;
static int CMISSProblemAdvectionDiffusionEquationType = 6;
static int CMISSProblemReactionDiffusionEquationType = 7;
static int CMISSProblemBiharmonicEquationType = 8;
static int CMISSProblemMonodomainEquationType = 1;
static int CMISSProblemBidomainEquationType = 2;
static int CMISSProblemLinearElasticModalType = 1;
static int CMISSProblemGalerkinProjectionType = 9;
static int CMISSProblemFiniteElasticityDarcyType = 1;
static int CMISSProblemFiniteElasticityStokesType = 2;
static int CMISSProblemFiniteElasticityNavierStokesType = 3;
static int CMISSProblemDiffusionDiffusionType = 4;
static int CMISSProblemDiffusionAdvectionDiffusionType = 5;
/*
 * >@}
 * > \addtogroup OPENCMISS_ProblemSubTypes OPENCMISS::Problem::Subtypes
 * > \brief Problem Subtypes.
 * > \see OPENCMISS::Problem,OPENCMISS
 * >@{
 */
static int CMISSProblemNoSubtype = 0;
static int CMISSProblemStaticStokesSubtype = 1;
static int CMISSProblemLaplaceStokesSubtype = 2;
static int CMISSProblemTransientStokesSubtype = 3;
static int CMISSProblemALEStokesSubtype = 5;
static int CMISSProblemOptimisedStokesSubtype = 4;
static int CMISSProblemStaticNavierStokesSubtype = 1;
static int CMISSProblemLaplaceNavierStokesSubtype = 2;
static int CMISSProblemTransientNavierStokesSubtype = 3;
static int CMISSProblemALENavierStokesSubtype = 5;
static int CMISSProblemOptimisedNavierStokesSubtype = 4;
static int CMISSProblemStandardDarcySubtype = 1;
static int CMISSProblemQuasistaticDarcySubtype = 2;
static int CMISSProblemALEDarcySubtype = 3;
static int CMISSProblemStandardLaplaceSubtype = 1;
static int CMISSProblemGeneralisedLaplaceSubtype = 2;
static int CMISSProblemLinearSourcePoissonSubtype = 1;
static int CMISSProblemNonlinearSourcePoissonSubtype = 2;
static int CMISSProblemNoSourceHelmholtzSubtype = 1;
static int CMISSProblemNoSourceDiffusionSubtype = 1;
static int CMISSProblemLinearSourceDiffusionSubtype = 2;
static int CMISSProblemNonlinearSourceDiffusionSubtype = 3;
static int CMISSProblemNoSourceAdvectionDiffusionSubtype = 1;
static int CMISSProblemLinearSourceAdvectionDiffusionSubtype = 2;
static int CMISSProblemNonlinearSourceAdvectionDiffusionSubtype = 3;
static int CMISSProblemNoSourceStaticAdvecDiffSubtype = 4;
static int CMISSProblemLinearSourceStaticAdvecDiffSubtype = 5;
static int CMISSProblemNonlinearSourceStaticAdvecDiffSubtype = 6;
static int CMISSProblemStandardGalerkinProjectionSubtype = 1;
static int CMISSProblemGeneralisedGalerkinProjectionSubtype = 2;
static int CMISSProblemMatPropertiesGalerkinProjectionSubtype = 3;
static int CMISSProblemStandardElasticityDarcySubtype = 101;
static int CMISSProblemCoupledSourceDiffusionDiffusionSubtype = 111;
static int CMISSProblemCoupledSourceDiffusionAdvecDiffusionSubtype = 121;
/*
 * >@}
 * > \addtogroup OPENCMISS_ProblemControlLoopTypes OPENCMISS::Problem::ControlLoopTypes
 * > \brief Problem control loop type parameters
 * > \see OPENCMISS::Problem,OPENCMISS
 * >@{
 */
static int CMISSProblemControlSimpleType = 1;
static int CMISSProblemControlFixedLoopType = 2;
static int CMISSProblemControlTimeLoopType = 3;
static int CMISSProblemControlWhileLoopType = 4;
/*
 * >@}
 * >@}
 */

/*
 * ==================================================================================================================================
 *
 *  PROBLEM_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * ==================================================================================================================================
 *
 *  REGION_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */


/*
 * ==================================================================================================================================
 *
 *  SOLVER_ROUTINES
 *
 * ==================================================================================================================================
 */

/* Module parameters */

/*
 * > \addtogroup OPENCMISS_SolverConstants OPENCMISS::Solver::Constants
 * > \brief Solver constants.
 * >@{
 * > \addtogroup OPENCMISS_SolverTypes OPENCMISS::Solver::SolverTypes
 * > \brief The types of solver
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverLinearType = 1;
static int CMISSSolverNonlinearType = 2;
static int CMISSSolverDynamicType = 3;
static int CMISSSolverDAEType = 4;
static int CMISSSolverEigenproblemType = 5;
static int CMISSSolverOptimiserType = 6;
/* >@}
 * > \addtogroup OPENCMISS_LinearSolverTypes OPENCMISS::Solver::LinearSolverTypes
 * > \brief The types of linear solvers.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverLinearDirectSolveType = 1;
static int CMISSSolverLinearIterativeSolveType = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_DirectLinearSolverTypes OPENCMISS::Solver::DirectLinearSolverTypes
 * > \brief The types of direct linear solvers. \todo Move libraries to a more appropriate place.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDirectLU = 1;
static int CMISSSolverDirectCholesky = 2;
static int CMISSSolverDirectSVD = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_IterativeLinearSolverTypes OPENCMISS::Solver::IterativeLinearSolverTypes
 * > \brief The types of iterative linear solvers.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverIterativeRichardson = 1;
static int CMISSSolverIterativeChebychev = 2;
static int CMISSSolverIterativeConjugateGradient = 3;
static int CMISSSolverIterativeBiconjugateGradient = 4;
static int CMISSSolverIterativeGMRES = 5;
static int CMISSSolverIterativeBiCGSTAB = 6;
static int CMISSSolverConjgradSquared = 7;
/*
 * >@}
 * > \addtogroup OPENCMISS_IterativePreconditionerTypes OPENCMISS::Solver::IterativePreconditionerTypes
 * > \brief The types of iterative preconditioners.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverIterativeNoPreconditioner = 0;
static int CMISSSolverIterativeJacobiPreconditioner = 1;
static int CMISSSolverIterativeBlockJacobiPreconditioner = 2;
static int CMISSSolverIterativeSORPreconditioner = 3;
static int CMISSSolverIterativeIncompleteCholeskyPreconditioner = 4;
static int CMISSSolverIterativeIncompleteLUPreconditioner = 5;
static int CMISSSolverIterativeAdditiveSchwarzPreconditioner = 6;
/*
 * >@}
 * > \addtogroup OPENCMISS_NonlinearSolverTypes OPENCMISS::Solver::NonlinearSolverTypes
 * > \brief The types of nonlinear solvers.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverNonlinearNewton = 1;
static int CMISSSolverNonlinearBFGSInverse = 2;
static int CMISSSolverNonlinearSQP = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_NewtonSolverTypes OPENCMISS::Solver::NewtonSolverTypes
 * > \brief The types of nonlinear Newton solvers.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverNewtonLinesearch = 1;
static int CMISSSolverNewtonTrustregion = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_NewtonLineSearchTypes OPENCMISS::Solver::NewtonLineSearchTypes
 * > \brief The types line search techniques for Newton line search nonlinear solvers.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverNewtonLinesearchNoNorms = 1;
static int CMISSSolverNewtonLinesearchNone = 2;
static int CMISSSolverNewtonLinesearchQuadratic = 3;
static int CMISSSolverNewtonLinesearchCubic = 4;
/*
 * >@}
 * > \addtogroup OPENCMISS_JacobianCalculationTypes OPENCMISS::Solver::JacobianCalculationTypes
 * > \brief The Jacobian calculation types for a nonlinear solver.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverNewtonJacobianNotCalculated = 1;
static int CMISSSolverNewtonJacobianAnalyticCalculated = 2;
static int CMISSSolverNewtonJacobianFDCalculated = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_DynamicOrderTypes OPENCMISS::Solver::DynamicOrderTypes
 * > \brief The order types for a dynamic solver.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDynamicFirstOrder = 1;
static int CMISSSolverDynamicSecondOrder = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_DynamicLinearityTypes OPENCMISS::Solver::DynamicLinearityTypes
 * > \brief The time linearity types for a dynamic solver.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDynamicLinear = 1;
static int CMISSSolverDynamicNonLinear = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_DynamicDegreeTypes OPENCMISS::Solver::DynamicDegreeTypes
 * > \brief The time interpolation polynomial degree types for a dynamic solver.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDynamicFirstDegree = 1;
static int CMISSSolverDynamicSecondDegree = 2;
static int CMISSSolverDynamicThirdDegree = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_DynamicSchemeTypes OPENCMISS::Solver::DynamicSchemeTypes
 * > \brief The types of dynamic solver scheme.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDynamicEulerScheme = 1;
static int CMISSSolverDynamicBackwardEulerScheme = 2;
static int CMISSSolverDynamicCrankNicholsonScheme = 3;
static int CMISSSolverDynamicGalerkinScheme = 4;
static int CMISSSolverDynamicZlamalScheme = 5;
static int CMISSSolverDynamicSecondDegreeGearScheme = 6;
static int CMISSSolverDynamicSecondDegreeLiniger1Scheme = 7;
static int CMISSSolverDynamicSecondDegreeLiniger2Scheme = 8;
static int CMISSSolverDynamicNewmark1Scheme = 9;
static int CMISSSolverDynamicNewmark2Scheme = 10;
static int CMISSSolverDynamicNewmark3Scheme = 11;
static int CMISSSolverDynamicThirdDegreeGearScheme = 12;
static int CMISSSolverDynamicThirdDegreeLiniger1Scheme = 13;
static int CMISSSolverDynamicThirdDegreeLiniger2Scheme = 14;
static int CMISSSolverDynamicHouboltScheme = 15;
static int CMISSSolverDynamicWilsonScheme = 16;
static int CMISSSolverDynamicBossakNewmark1Scheme = 17;
static int CMISSSolverDynamicBossakNewmark2Scheme = 18;
static int CMISSSolverDynamicHilbertHughesTaylor1Scheme = 19;
static int CMISSSolverDynamicHilbertHughesTaylor2Scheme = 20;
static int CMISSSolverDynamicUserDefinedScheme = 21;
/*
 * >@}
 * > \addtogroup OPENCMISS_DAETypes OPENCMISS::Solver::DAETypes
 * > \brief The type of differential-algebraic equation.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDAEDifferentialOnly = 0;
static int CMISSSolverDAEIndex1 = 1;
static int CMISSSolverDAEIndex2 = 2;
static int CMISSSolverDAEIndex3 = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_DAESolverTypes OPENCMISS::Solver::DAESolverTypes
 * > \brief The differential-algebraic equation solver types for a differential-algebraic equation solver.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDAEEuler = 1;
static int CMISSSolverDAECrankNicholson = 2;
static int CMISSSolverDAERungeKutta = 3;
static int CMISSSolverDAEAdamasMoulton = 4;
static int CMISSSolverDAEBDF = 5;
static int CMISSSolverRushLarson = 6;
/*
 * >@}
 * > \addtogroup OPENCMISS_EulerDAESolverTypes OPENCMISS::Solver::EulerDAESolverTypes
 * > \brief The Euler solver types for a differential-algebriac equation solver.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverDAEEulerForward = 1;
static int CMISSSolverDAEEulerBackward = 2;
static int CMISSSolverDAEEulerImproved = 3;
/*
 * >@}
 * > \addtogroup OPENCMISS_SolutionInitialiseTypes OPENCMISS::Solver::SolutionInitialiseTypes
 * > \brief The types of solution initialisation.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverSolutionInitialiseZero = 0;
static int CMISSSolverSolutionInitialiseCurrentField = 1;
static int CMISSSolverSolutionInitialiseNoChange = 2;
/*
 * >@}
 * > \addtogroup OPENCMISS_SolverOutputTypes OPENCMISS::Solver::OutputTypes
 * > \brief The types of output.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverNoOutput = 0;
static int CMISSSolverProgressOutput = 1;
static int CMISSSolverTimingOutput = 2;
static int CMISSSolverSolverOutput = 3;
static int CMISSSolverSolverMatrixOutput = 4;
/*
 * >@}
 * > \addtogroup OPENCMISS_SolverEquationsSparsityTypes OPENCMISS::SolverEquations::SparsityTypes
 * > \brief The types of sparse solver equations matrices.
 * > \see OPENCMISS::Solver::Constants,OPENCMISS
 * >@{
 */
static int CMISSSolverEquationsSparseMatrices = 1;
static int CMISSSolverEquationsFullMatrices = 2;
/*
 * >@}
 * >@}
 */

















/*
 *=================================================================================================================================
 *
 * CMISS
 *
 *==================================================================================================================================
 */

/* 
 * Struct defs
 */

struct CMISSBasisType_;
struct CMISSBoundaryConditionsType_;
struct CMISSCellMLType_;
struct CMISSControlLoopType_;
struct CMISSCoordinateSystemType_;
struct CMISSDecompositionType_;
struct CMISSEquationsType_;
struct CMISSEquationsSetType_;
struct CMISSFieldType_;
struct CMISSFieldsType_;
struct CMISSGeneratedMeshType_;
struct CMISSHistoryType_;
struct CMISSMeshType_;
struct CMISSMeshElementsType_;
struct CMISSNodesType_;
struct CMISSProblemType_;
struct CMISSQuadratureType_;
struct CMISSRegionType_;
struct CMISSSolverType_;
struct CMISSSolverEquationsType_;

/* 
 * Type defs
 */



typedef int CMISSError;
typedef struct CMISSBasisType_ *CMISSBasisType;
typedef struct CMISSBoundaryConditionsType_ *CMISSBoundaryConditionsType;
typedef struct CMISSCellMLType_ *CMISSCellMLType;
typedef struct CMISSControlLoopType_ *CMISSControlLoopType;
typedef struct CMISSCoordinateSystemType_ *CMISSCoordinateSystemType;
typedef struct CMISSDecompositionType_ *CMISSDecompositionType;
typedef struct CMISSEquationsType_ *CMISSEquationsType;
typedef struct CMISSEquationsSetType_ *CMISSEquationsSetType;
typedef struct CMISSFieldType_ *CMISSFieldType;
typedef struct CMISSFieldsType_ *CMISSFieldsType;
typedef struct CMISSGeneratedMeshType_ *CMISSGeneratedMeshType;
typedef struct CMISSHistoryType_ *CMISSHistoryType;
typedef struct CMISSMeshType_ *CMISSMeshType;
typedef struct CMISSMeshElementsType_ *CMISSMeshElementsType;
typedef struct CMISSNodesType_ *CMISSNodesType;
typedef struct CMISSProblemType_ *CMISSProblemType;
typedef struct CMISSQuadratureType_ *CMISSQuadratureType;
typedef struct CMISSRegionType_ *CMISSRegionType;
typedef struct CMISSSolverType_ *CMISSSolverType;
typedef struct CMISSSolverEquationsType_ *CMISSSolverEquationsType;

/* 
 * Protypes
 */

/*
 *=================================================================================================================================
 *
 * Types
 *
 *==================================================================================================================================
 */

CMISSError CMISSBasisTypeFinalise(CMISSBasisType *Basis);

CMISSError CMISSBasisTypeInitialise(CMISSBasisType *Basis);

CMISSError CMISSBoundaryConditionsTypeFinalise(CMISSBoundaryConditionsType *BoundaryConditions);

CMISSError CMISSBoundaryConditionsTypeInitialise(CMISSBoundaryConditionsType *BoundaryConditions);

CMISSError CMISSCellMLTypeFinalise(CMISSCellMLType *CellML);

CMISSError CMISSCellMLTypeInitialise(CMISSCellMLType *CellML);

CMISSError CMISSControlLoopTypeFinalise(CMISSControlLoopType *ControlLoop);

CMISSError CMISSControlLoopTypeInitialise(CMISSControlLoopType *ControlLoop);

CMISSError CMISSCoordinateSystemTypeFinalise(CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSCoordinateSystemTypeInitialise(CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSDecompositionTypeFinalise(CMISSDecompositionType *Decomposition);

CMISSError CMISSDecompositionTypeInitialise(CMISSDecompositionType *Decomposition);

CMISSError CMISSEquationsTypeFinalise(CMISSEquationsType *Equations);

CMISSError CMISSEquationsTypeInitialise(CMISSEquationsType *Equations);

CMISSError CMISSEquationsSetTypeFinalise(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetTypeInitialise(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSFieldTypeFinalise(CMISSFieldType *Field);

CMISSError CMISSFieldTypeInitialise(CMISSFieldType *Field);

CMISSError CMISSFieldsTypeCreate(const CMISSRegionType Region,
		CMISSFieldsType *Fields);

CMISSError CMISSFieldsTypeFinalise(CMISSFieldsType *Fields);

CMISSError CMISSFieldsTypeInitialise(CMISSFieldsType *Fields);

CMISSError CMISSGeneratedMeshTypeFinalise(CMISSGeneratedMeshType *GeneratedMesh);

CMISSError CMISSGeneratedMeshTypeInitialise(CMISSGeneratedMeshType *GeneratedMesh);

CMISSError CMISSHistoryTypeFinalise(CMISSHistoryType *History);

CMISSError CMISSHistoryTypeInitialise(CMISSHistoryType *History);

CMISSError CMISSMeshTypeFinalise(CMISSMeshType *Mesh);

CMISSError CMISSMeshTypeInitialise(CMISSMeshType *Mesh);

CMISSError CMISSMeshElementsTypeFinalise(CMISSMeshElementsType *MeshElements);

CMISSError CMISSMeshElementsTypeInitialise(CMISSMeshElementsType *MeshElements);

CMISSError CMISSNodesTypeFinalise(CMISSNodesType *Nodes);

CMISSError CMISSNodesTypeInitialise(CMISSNodesType *Nodes);

CMISSError CMISSProblemTypeFinalise(CMISSProblemType *Problem);

CMISSError CMISSProblemTypeInitialise(CMISSProblemType *Problem);

CMISSError CMISSQuadratureTypeFinalise(CMISSQuadratureType *Quadrature);

CMISSError CMISSQuadratureTypeInitialise(CMISSQuadratureType *Quadrature);

CMISSError CMISSRegionTypeFinalise(CMISSRegionType *Region);

CMISSError CMISSRegionTypeInitialise(CMISSRegionType *Err);

CMISSError CMISSSolverTypeFinalise(CMISSSolverType *Solver);

CMISSError CMISSSolverTypeInitialise(CMISSSolverType *Solver);

CMISSError CMISSSolverEquationsTypeFinalise(CMISSSolverEquationsType *SolverEquations);

CMISSError CMISSSolverEquationsTypeInitialise(CMISSSolverEquationsType *SolverEquations);

/*
 *=================================================================================================================================
 *
 * ANALYTIC_ANALYSIS_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSAnalyticAnalysisOutputNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int FileNameSize,
		const char *FileName);

CMISSError CMISSAnalyticAnalysisOutput(const CMISSFieldType Field,
		const int FileNameSize,
		const char *FileName);

CMISSError CMISSAnalyticAnalysisAbsoluteErrorGetNodeNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int DerivativeNumber,
		const int NodeNumber,
		const int ComponentNumber,
		const int VariableNumber,
		double *Value);

CMISSError CMISSAnalyticAnalysisAbsoluteErrorGetNode(const CMISSFieldType Field,
		const int DerivativeNumber,
		const int NodeNumber,
		const int ComponentNumber,
		const int VariableNumber,
		double *Value);

CMISSError CMISSAnalyticAnalysisPercentageErrorGetNodeNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int DerivativeNumber,
		const int NodeNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisPercentageErrorGetNode(const CMISSFieldType Field,
		const int DerivativeNumber,
		const int NodeNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRelativeErrorGetNodeNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int DerivativeNumber,
		const int NodeNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRelativeErrorGetNode(const CMISSFieldType Field,
		const int DerivativeNumber,
		const int NodeNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisAbsoluteErrorGetElementNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ElementNumber,
		const int ComponentNumber,
		const int VariableType,
		doulbe *Value);

CMISSError CMISSAnalyticAnalysisAbsoluteErrorGetElement(const CMISSFieldType Field,
		const int ElementNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisPercentageErrorGetElementNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ElementNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisPercentageErrorGetElement(const CMISSFieldType Field,
		const int ElementNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRelativeErrorGetElementNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ElementNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRelativeErrorGetElement(const CMISSFieldType Field,
		const int ElementNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisAbsoluteErrorGetConstantNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisAbsoluteErrorGetConstant(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisPercentageErrorGetConstantNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableTYpe,
		double *Value);

CMISSError CMISSAnalyticAnalysisPercentageErrorGetConstant(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRelativeErrorGetConstantNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRelativeErrorGetConstant(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *Value);

CMISSError CMISSAnalyticAnalysisRmsErrorGetNodeNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		const int ErrorType,
		double *LocalValue,
		double *LocalGhostValue)
		double *GlobalValue);

CMISSError CMISSAnalyticAnalysisRmsErrorGetNode(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		const int ErrorType,
		double *LocalValue,
		double *LocalGhostValue,
		double *GlobalValue);

CMISSError CMISSAnalyticAnalysisRmsErrorGetElementNumber(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		const int ErrorType,
		double *LocalValue,
		double *LocalGhostValue)
		double *GlobalValue);

CMISSError CMISSAnalyticAnalysisRmsErrorGetElement(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		const int ErrorType,
		double *LocalValue,
		double *LocalGhostValue,
		double *GlobalValue);

CMISSError CMISSAnalyticAnalysisIntegralNumericalValueGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralNumericalValueGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralAnalyticValueGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralAnalyticValueGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralPercentageErrorGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralPercentageErrorGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralAbsoluteErrorGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralAbsoluteErrorGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralRelativeErrorGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralRelativeErrorGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralNidNumericalValueGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralNidNumericalValueGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralNidErrorGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

CMISSError CMISSAnalyticAnalysisIntegralNidErrorGet(const CMISSFieldType Field,
		const int ComponentNumber,
		const int VariableType,
		double *IntegralValue,
		double *GhostIntegralValue);

/*
 *=================================================================================================================================
 *
 * BASE_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSDiagnosticsSetOff();

CMISSError CMISSDiagnosticsSetOn(const int DiagType,
		const int LevelListSize[1],
		const int *LevelList,
		const int DiagFilenameSize,
		const char *DiagFilename,
		const int RoutineListSize[1],
		const char *RoutineList);

CMISSError CMISSOutputSetOff();

CMISSError CMISSOutputSetOn(const int EchoFilenameSize,
		const char *EchoFilename);

CMISSError CMISSTimingSetOff();

CMISSError CMISSTimingSetOn(const int TimingType,
		const int TimingSummaryFlag,
		const char *TimingFilename,
		const int RoutineListSize[1],
		const char *RoutineList);

CMISSError CMISSTimingSummaryOutput();

/*
 *==================================================================================================================================
 *
 * BASIS_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSBasisCollapsedXiGetNum(const int UserNumber,
		int *CollapsedXiSize,
		int *CollapsedXi);

CMISSError CMISSBasisCollapsedXiGet(const CMISSBasisType Basis,
		int *CollapsedXiSize,
		int *CollapsedXi);

CMISSError CMISSBasisCollapsedXiSetNum(const int UserNumber,
		const int CollapsedXiSize[1],
		const int *CollapsedXi);

CMISSError CMISSBasisCollapsedXiSet(CMISSBasisType *Basis,
		const int CollapsedXiSize[1],
		const int *CollapsedXi);

CMISSError CMISSBasisCreateFinishNum(const int UserNumber);

CMISSError CMISSBasisCreateFinish(CMISSBasisType *Basis);

CMISSError CMISSBasisCreateStartNum(const int UserNumber);

CMISSError CMISSBasisCreateStart(const int UserNumber,
		CMISSBasisType *Basis);

CMISSError CMISSBasisDestroyNum(const int UserNumber);

CMISSError CMISSBasisDestroy(CMISSBasisType *Basis);

CMISSError CMISSBasisInterpolationXiGetNum(const int UserNumber,
		int *InterpolationXiSize,
		int *InterpolationXi);

CMISSError CMISSBasisInterpolationXiGet(const CMISSBasisType Basis,
		int *InterpolationXiSize,
		int *InterpolationXi);

CMISSError CMISSBasisInterpolationXiSetNum(const int UserNumber,
		const int InterpolationXiSize[1],
		const int *InterpolationXi);

CMISSError CMISSBasisInterpolationXiSet(const CMISSBasisType Basis,
		const int InterpolationXiSize[1],
		const int *InterpolationXi);

CMISSError CMISSBasisNumberOfLocalNodesGetNum(const int UserNumber,
		int *NumberOfLocalNodes);

CMISSError CMISSBasisNumberOfLocalNodesGet(const CMISSBasisType Basis,
		int *NumberOfNodes);

CMISSError CMISSBasisNumberOfXiGetNum(const int UserNumber,
		int *NumberOfXi);

CMISSError CMISSBasisNumberOfXiGet(const CMISSBasisType Basis,
		int *NumberOfXi);

CMISSError CMISSBasisNumberOfXiSetNum(const int UserNumber,
		const int NumberOfXi);

CMISSError CMISSBasisNumberOfXiSet(CMISSBasisType *Basis,
		const int NumberOfXi);

CMISSError CMISSBasisQuadratureNumberOfGaussXiGetNum(const int UserNumber,
		int *NumberOfGaussXiSize,
		int *NumberOfGaussXi);

CMISSError CMISSBasisQuadratureNumberOfGaussXiGet(const CMISSBasisType Basis,
		int *NumberOfGaussXiSize,
		int *NumberOfGaussXi);

CMISSError CMISSBasisQuadratureNumberOfGaussXiSetNum(const int UserNumber,
		const int NumberOfGaussXiSize[1],
		const int *NumberOfGaussXi);

CMISSError CMISSBasisQuadratureNumberOfGaussXiSet(CMISSBasisType *Basis,
		const int NumberOfGaussXiSize[1],
		const int *NumberOfGaussXi);

CMISSError CMISSBasisQuadratureOrderGetNum(const int UserNumber,
		int *QuadratureOrder);

CMISSError CMISSBasisQuadratureOrderGet(const CMISSBasisType Basis,
		int *QuadratureOrder);

CMISSError CMISSBasisQuadratureOrderSetNum(const int UserNumber,
		const int QuadratureOrder);

CMISSError CMISSBasisQuadratureOrderSet(CMISSBasisType *Basis,
		const int QuadratureOrder);

CMISSError CMISSBasisQuadratureTypeGetNum(const int UserNumber,
		int *QuadratureType);

CMISSError CMISSBasisQuadratureTypeGet(const CMISSBasisType Basis,
		int *QuadratureType);

CMISSError CMISSBasisQuadratureTypeSetNum(const int UserNumber,
		const int QuadratureType);

CMISSError CMISSBasisQuadratureTypeSet(CMISSBasisType *Basis,
		const int QuadratureType);

CMISSError CMISSBasisTypeGetNum(const int UserNumber,
		int *BasisType);

CMISSError CMISSBasisTypeGet(const CMISSBasisType Basis,
		int *BasisType);

CMISSError CMISSBasisTypeSetNum(const int UserNumber,
		const int BasisType);

CMISSError CMISSBasisTypeSet(CMISSBasisType *Basis,
		const int BasisType);

/*
 *==================================================================================================================================
 *
 * BOUNDARY_CONDITIONS_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSBoundaryConditionsDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSBoundaryConditionsDestroy(CMISSBoundaryConditionsType *BoundaryConditions);

CMISSError CMISSBoundaryConditionsAddConstantNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int VariableType,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsAddConstant(const CMISSBoundaryConditionsType BoundaryConditions,
		const int VariableType,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsSetConstantNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int VariableType,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsSetConstant(const CMISSBoundaryConditionsType BoundaryConditions,
		const int VariableType,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsAddElementNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int VariableType,
		const int ElementUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsAddElement(const CMISSBoundaryConditionsType BoundaryConditions,
		const int VariableType,
		const int ElementUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsSetElementNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int VariableType,
		const int ElementUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsSetElement(const CMISSBoundaryConditionsType BoundaryConditions,
		const int VariableType,
		const int ElementUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsAddNodeNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int VariableType,
		const int DerivativeNumber,
		const int NodeUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsAddNode(const CMISSBoundaryConditionsType BoundaryConditions,
		const int VariableType,
		const int DerivativeNumber,
		const int NodeUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsSetNodeNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int VariableType,
		const int DerivativeNumber,
		const int NodeUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSBoundaryConditionsSetNode(const CMISSBoundaryConditionsType BoundaryConditions,
		const int VariableType,
		const int DerivativeNumber,
		const int NodeUserNumber,
		const int ComponentNumber,
		const int Condition,
		const double Value);

CMISSError CMISSEquationsSetBoundaryConditionsGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		CMISSBoundaryConditionsType *BoundaryConditions);

CMISSError CMISSEquationsSetBoundaryConditionsGet(const CMISSEquationsSetType EquationsSet,
		CMISSBoundaryConditionsType *BoundaryConditions);




/*
 *=================================================================================================================================
 *
 * CMISS
 *
 *==================================================================================================================================
 */

CMISSError CMISSErrorHandlingModeGet(int *ErrorHandlingMode);

CMISSError CMISSErrorHandlingModeSet(const int ErrorHandlingMode);

CMISSError CMISSFinalise();

CMISSError CMISSInitialiseNum(int *WorldCoordinateSystemUserNumber,
			      int *WorldRegionUserNumber);

CMISSError CMISSInitialise(CMISSCoordinateSystemType *WorldCoordinateSystem,
			   CMISSRegionType *WorldRegion);

/*
 *==================================================================================================================================
 *
 * CMISS_CELLML
 *
 *==================================================================================================================================
 */

CMISSError CMISSCellMLCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLCreateFinish(CMISSCellMLType *CellML)

CMISSError CMISSCellMLCreateStartNum(const int CellMLUserNumber,
		const int RegionUserNumber,
		const int FieldUserNumber);

CMISSError CMISSCellMLCreateStart(const int CellMLUserNumber,
		CMISSFieldType *Field,
		CMISSCellMLType *CellML);

CMISSError CMISSCellMLDestroyNum(const int CellMLUserNumber);

CMISSError CMISSCellMLDestroy(CMISSCellMLType *CellML);

CMISSError CMISSCellMLModelsCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLModelsCreateFinish(CMISSCellMLType *CellML);

CMISSError CMISSCellMLModelsCreateStartNum(const int CellMLUserNumber);

CMISSError CMISSCellMLModelsCreateStart(CMISSCellMLType *CellML);

CMISSError CMISSCellMLModelImportNum(const int CellMLModelUserNumber,
		const int CellMLUserNumber,
		const int URISize,
		char *URI);

CMISSError CMISSCellMLModelImport(const int CellMLModelUserNumber,
		CMISSCellMLType *CellML,
		const int URISize,
		char *URI);

CMISSError CMISSCellMLModelsFieldCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLModelsFieldCreateFinish(CMISSCellMLType *CellML);

CMISSError CMISSCellMLModelsFieldCreateStartNum(const int CellMLModelsFieldUserNumber,
		const int CellMLUserNumber);

CMISSError CMISSCellMLModelsFieldCreateStart(const int CellMLModelsFieldUserNumber,
		CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLModelsFieldGetNum(const int CellMLUserNumber,
		int *CellMLModelsFieldUserNumber);

CMISSError CMISSCellMLModelsFieldGet(CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLStateFieldCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLStateFieldCreateFinish(CMISSCellMLType *CellML);

CMISSError CMISSCellMLStateFieldCreateStartNum(const int CellMLStateFieldUserNumber,
		const int CellMLUserNumber);

CMISSError CMISSCellMLStateFieldCreateStart(const int CellMLStateFieldUserNumber,
		CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLStateFieldGetNum(const int CellMLUserNumber,
		int *CellMLStateFieldUserNumber);

CMISSError CMISSCellMLStateFieldGet(CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLFieldComponentGetNum(const int CellMLUserNumber,
		const int CellMLFieldType,
		const int URISize,
		const char *URI,
		int *FieldComponent);

CMISSError CMISSCellMLFieldComponentGet(CMISSCellMLType *CellML,
		const int CellMLFieldType,
		const int URISize,
		const char *URI,
		int *FieldComponent);

CMISSError CMISSCellMLIntermediateFieldAddNum(const int CellMLUserNumber,
		const int CellMLModelUserNumber,
		const int URISize,
		const char *URI);

CMISSError CMISSCellMLIntermediateFieldAdd(CMISSCellMLType *CellML,
		const int CellMLModelUserNumber,
		const int URISize,
		const char *URI);

CMISSError CMISSCellMLIntermediateFieldCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLIntermediateFieldCreateFinish(CMISSCellMLType *CellML);

CMISSError CMISSCellMLIntermediateFieldCreateStartNum(const int CellMLIntermediateFieldUserNumber,
		const int CellMLUserNumber);

CMISSError CMISSCellMLIntermediateFieldCreateStart(const int CellMLIntermediateFieldUserNumber,
		CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLIntermediateFieldGetNum(const int CellMLUserNumber,
		int *CellMLIntermediateFieldUserNumber);

CMISSError CMISSCellMLIntermediateFieldGet(CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLParameterAddNum(const int CellMLUserNumber,
	const int CellMLModelUserNumber,
	const int URISize,
	const char *URI);

CMISSError CMISSCellMLParameterAdd(CMISSCellMLType *CellML,
	const int CellMLModelUserNumber,
	const int URISize,
	const char *URI);

CMISSError CMISSCellMLParametersCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLParametersCreateFinish(CMISSCellMLType *CellML);

CMISSError CMISSCellMLParametersCreateStartNum(const int CellMLUserNumber);

CMISSError CMISSCellMLParametersCreateStart(CMISSCellMLType *CellML);

CMISSError CMISSCellMLParametersFieldCreateFinishNum(const int CellMLUserNumber);

CMISSError CMISSCellMLParametersFieldCreateFinish(CMISSCellMLType *CellML);

CMISSError CMISSCellMLParametersFieldCreateStartNum(const int CellMLParametersFieldUserNumber,
		const int CellMLUserNumber);

CMISSError CMISSCellMLParametersFieldCreateStart(const int CellMLParametersFieldUserNumber,
		CMISSCellMLType *CellML,
		CMISSFieldType *Field);

CMISSError CMISSCellMLParametersFieldGetNum(const int CellMLUserNumber,
		int *CellMLParametersFieldUserNumber);

CMISSError CMISSCellMLParametersFieldGet(CMISSCellMLType *CellML.
		CMISSFieldType *Field);

CMISSError CMISSCellMLGenerateNum(const int CellMLUserNumber);

CMISSError CMISSCellMLGenerate(CMISSCellMLType *CellML);

/*
 *==================================================================================================================================
 *
 * COMP_ENVIRONMENT
 *
 *==================================================================================================================================
 */

CMISSError CMISSComputationalNodeNumberGet(int *NodeNumber);

CMISSError CMISSComputationalNumberOfNodesGet(int *NumberOfNodes);

/*
 *==================================================================================================================================
 *
 * CONTROL_LOOP_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSControlLoopCurrentTimesGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		double *CurrentTime,
		double *TimeIncrement);

CMISSError CMISSControlLoopCurrentTimesGet(const CMISSControlLoopType ControlLoop,
		double *CurrentTime,
		double *TimeIncrement);

CMISSError CMISSControlLoopDestroyNum(const int ProblemUserNumber,
		const int ControlLoopIdentifitersSize[1],
		const int *ControlLoopIdentifiers);

CMISSError CMISSControlLoopDestroy(CMISSControlLoopType *ControlLoop);

CMISSError CMISSControlLoopGetNum(const int ProblemUserNumber,
		const int ControlLoopRootIdentifiersSize[1],
		const int *ControlLoopRootIdentifiers,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		CMISSControlLoopType *ControlLoop);

CMISSError CMISSControlLoopGet(const CMISSControlLoopType ControlLoopRoot,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		CMISSControlLoopType *ControlLoop);

CMISSError CMISSControlLoopIterationsSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int StartIteration,
		const int StopIteration,
		const int IterationIncrement);

CMISSError CMISSControlLoopIterationsSet(CMISSControlLoopType *ControlLoop,
		const int StartIteration,
		const int StopIteration,
		const int IterationIncrement);

CMISSError CMISSControlLoopMaximumIterationsSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int MaximumIterations);

CMISSError CMISSControlLoopMaximumIterationsSet(CMISSControlLoopType *ControlLoop,
		const int MaximumIterations);

CMISSError CMISSControlLoopNumberOfSubLoopsGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		int *NumberOfSubLoops);

CMISSError CMISSControlLoopNumberOfSubLoopsGet(const CMISSControlLoopType ControlLoop,
		int *NumberOfSubLoops);

CMISSError CMISSControlLoopNumberOfSubLoopsSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int NumberOfSubLoops);

CMISSError CMISSControlLoopNumberOfSubLoopsSet(CMISSControlLoopType *ControlLoop,
		const int NumberOfSubLoops);

CMISSError CMISSControlLoopTimeOutputSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int OutputFrequency);

CMISSError CMISSControlLoopTimeOutputSet(CMISSControlLoopType *ControlLoop,
		const int OutputFrequency);

CMISSError CMISSControlLoopTimeInputSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int InputOption);

CMISSError CMISSControlLoopTimeInputSet(CMISSControlLoopType *ControlLoop,
		const int InputOption);

CMISSError CMISSControlLoopTimesGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		double *StartTime,
		double *StopTime,
		double *TimeIncrement,
		double *CurrentTime);

CMISSError CMISSControlLoopTimesGet(const CMISSControlLoopType ControlLoop,
		double *StartTime,
		double *StopTime,
		double *TimeIncrement,
		double *CurrentTime);

CMISSError CMISSControlLoopTimesSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const double StartTime,
		const double StopTime,
		const double TimeIncrement);

CMISSError CMISSControlLoopTimesSet(CMISSControlLoopType *ControlLoop,
		const double StartTime,
		const double StopTime,
		const double TimeIncrement);

CMISSError CMISSControlLoopTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int LoopType);

CMISSError CMISSControlLoopTypeSet(CMISSControlLoopType *ControlLoop,
		const int LoopType);

/*
 * ==================================================================================================================================
 *
 * COORDINATE_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSCoordinateSystemCreateFinishNum(const int CoordinateSystemUserNumber);

CMISSError CMISSCoordinateSystemCreateFinish(CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSCoordinateSystemCreateStartNum(const int CoordinateSystemUserNumber);

CMISSError CMISSCoordinateSystemCreateStart(const int CoordinateSystemUserNumber,
		CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSCoordinateSystemDestroyNum(const int CoordinateSystemUserNumber);

CMISSError CMISSCoordinateSystemDestroy(CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSCoordinateSystemDimensionGetNum(const int CoordinateSystemUserNumber,
		int *CoordinateSystemDimension);

CMISSError CMISSCoordinateSystemDimensionGet(const CMISSCoordinateSystemType CoordinateSystem,
		int *CoordinateSystemDimension);

CMISSError CMISSCoordinateSystemDimensionSetNum(const int CoordinateSystemUserNumber,
		const int CoordinateSystemDimension);

CMISSError CMISSCoordinateSystemDimensionSet(CMISSCoordinateSystemType *CoordinateSystem,
		const int CoordinateSystemDimension);

CMISSError CMISSCoordinateSystemFocusGetNum(const int CoordinateSystemUserNumber,
		double *Focus);

CMISSError CMISSCoordinateSystemFocusGet(const CMISSCoordinateSystemType CoordinateSystem,
		double *Focus);

CMISSError CMISSCoordinateSystemFocusSetNum(const int CoordinateSystemUserNumber,
		const double Focus);

CMISSError CMISSCoordinateSystemFocusSet(CMISSCoordinateSystemType *CoordinateSystem,
		const double Focus);

CMISSError CMISSCoordinateSystemRadialInterpolationGetNum(const int CoordinateSystemUserNumber,
		int *RadialInterpolationType);

CMISSError CMISSCoordinateSystemRadialInterpolationGet(CMISSCoordinateSystemType *CoordinateSystem,
		int *RadialInterpolationType);

CMISSError CMISSCoordinateSystemRadialInterpolationSetNum(const int CoordinateSystemUserNumber,
		const int RadialInterpolationType);

CMISSError CMISSCoordinateSystemRadialInterpolationSet(CMISSCoordinateSystemType *CoordinateSystem,
		const int RadialInterpolationType);

CMISSError CMISSCoordinateSystemTypeGetNum(const int CoordinateSystemUserNumber,
		int *CoordinateSystemType);

CMISSError CMISSCoordinateSystemTypeGet(const CMISSCoordinateSystemType CoordinateSystem,
		int *CoordinateSystemType);

CMISSError CMISSCoordinateSystemTypeSetNum(const int CoordinateSystemUserNumber,
		const int CoordinateSystemType);

CMISSError CMISSCoordinateSystemTypeSet(CMISSCoordinateSystemType *CoordinateSystem,
		const int CoordinateSystemType);

CMISSError CMISSCoordinateSystemOriginGetNum(const int CoordinateSystemUserNumber,
		int *OriginSize[1],
		double *Origin);

CMISSError CMISSCoordinateSystemOriginGet(const CMISSCoordinateSystemType CoordinateSystem,
		int *OriginSize[1],
		double *Origin);

CMISSError CMISSCoordinateSystemOriginSetNum(const int CoordinateSystemUserNumber,
		const int OriginSize[1],
		const double *Origin);

CMISSError CMISSCoordinateSystemOriginSet(const CMISSCoordinateSystemType CoordinateSystem,
		const int OriginSize[1],
		const double *Origin);

CMISSError CMISSCoordinateSystemOrientationGetNum(const int CoordinateSystemUserNumber,
		int *OrientationSize[1],
		double *Orientation);

CMISSError CMISSCoordinateSystemOrientationGet(const CMISSCoordinateSystemType CoordinateSystem,
		int *OrientationSize[1],
		double *Orientation);

CMISSError CMISSCoordinateSystemOrientationSetNum(const int CoordinateSystemUserNumber,
		const int OrientationSize[1],
		const double *Orientation);

CMISSError CMISSCoordinateSystemOrientationSet(CMISSCoordinateSystemType CoordinateSystem,
		const int OrientationSize[1],
		const double *Orientation);

/*
 * ==================================================================================================================================
 *
 * EQUATIONS_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSEquationsDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsDestroy(CMISSEquationsType *Equations);

CMISSError CMISSEquationsLinearityTypeGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *LinearityType);

CMISSError CMISSEquationsLinearityTypeGet(const CMISSEquationsType Equations,
		int *LinearityType);

CMISSError CMISSEquationsLumpingTypeGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *LumpingType);

CMISSError CMISSEquationsLumpingTypeGet(const CMISSEquationsType Equations,
		int LumpingType);

CMISSError CMISSEquationsLumpingTypeSetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int LumpingType);

CMISSError CMISSEquationsLumpingTypeSet(CMISSEquationsType *Equations,
		const int LumpingType);

CMISSError CMISSEquationsOutputTypeGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *OutputType);

CMISSError CMISSEquationsOutputTypeGet(const CMISSEquationsType Equations,
		int *OutputType);

CMISSError CMISSEquationsOutputTypeSetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int OutputType);

CMISSError CMISSEquationsOutputTypeSet(CMISSEquationsType *Equations,
		const int OutputType);

CMISSError CMISSEquationsSparsityTypeGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *SparsityType);

CMISSError CMISSEquationsSparsityTypeGet(const CMISSEquationsType Equations,
		int *SparsityType);

CMISSError CMISSEquationsSparsityTypeSetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int SparsityType);

CMISSError CMISSEquationsSparsityTypeSet(CMISSEquationsType *Equations,
		const int SparsityType);

CMISSError CMISSEquationsTimeDependenceTypeGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *TimeDependenceType);

CMISSError CMISSEquationsTimeDependenceTypeGet(const CMISSEquationsType Equations,
		int *TimeDependenceType);

/*
 *==================================================================================================================================
 *
 * EQUATIONS_SET_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSEquationsSetAnalyticCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetAnalyticCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetAnalyticCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int AnalyticFunctionType,
		const int AnalyticFieldUserNumber);

CMISSError CMISSEquationsSetAnalyticCreateStart(CMISSEquationsSetType *EquationsSet,
		const int AnalyticFunctionType,
		const int AnalyticFieldUserNumber,
		CMISSFieldType *AnalyticField);

CMISSError CMISSEquationsSetAnalyticDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetAnalyticDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetBoundaryConditionsAnalyticNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetBoundaryConditionsAnalytic(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetBoundaryConditionsCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetBoundaryConditionsCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetBoundaryConditionsCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetBoundaryConditionsCreateStart(CMISSEquationsSetType *EquationsSet,
		CMISSBoundaryConditionsType *BoundaryConditions);

CMISSError CMISSEquationsSetBoundaryConditionsDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetBoundaryConditionsDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int GeomFibreFieldUserNumber);

CMISSError CMISSEquationsSetCreateStart(const int EquationsSetUserNumber,
		const CMISSRegionType Region,
		const CMISSFieldType GeomFibreField,
		CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetDependentCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetDependentCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetDependentCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int DependentFieldUserNumber);

CMISSError CMISSEquationsSetDependentCreateStart(CMISSEquationsSetType *EquationsSet,
		const int DependentFieldUserNumber,
		CMISSFieldType *DependentField);

CMISSError CMISSEquationsSetDependentDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetDependentDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetEquationsCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetEquationsCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetEquationsCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetEquationsCreateStart(CMISSEquationsSetType *EquationsSet,
		CMISSEquationsType *Equations);

CMISSError CMISSEquationsSetEquationsDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetEquationsDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetIndependentCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetIndependentCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetIndependentCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int IndependentFieldUserNumber);

CMISSError CMISSEquationsSetIndependentCreateStart(CMISSEquationsSetType *EquationsSet,
		const int IndependentFieldUserNumber,
		CMISSFieldType *IndependentField);

CMISSError CMISSEquationsSetIndependentDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetIndependentDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetMaterialsCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetMaterialsCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetMaterialsCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int MaterialsFieldUserNumber);

CMISSError CMISSEquationsSetMaterialsCreateStart(CMISSEquationsSetType *EquationsSet,
		const int MaterialsFieldUserNumber,
		CMISSFieldType *MaterialsField);

CMISSError CMISSEquationsSetMaterialsDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetMaterialsDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetSolutionMethodGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *SolutionMethod);

CMISSError CMISSEquationsSetSolutionMethodGet(CMISSEquationsSetType *EquationsSet,
		int *SolutionMethod);

CMISSError CMISSEquationsSetSolutionMethodSetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int SolutionMethod);

CMISSError CMISSEquationsSetSolutionMethodSet(CMISSEquationsSetType *EquationsSet,
		const int SolutionMethod);

CMISSError CMISSEquationsSetSourceCreateFinishNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetSourceCreateFinish(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetSourceCreateStartNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int SourceFieldUserNumber);

CMISSError CMISSEquationsSetSourceCreateStart(CMISSEquationsSetType *EquationsSet,
		const int SourceFieldUserNumber,
		CMISSFieldType *SourceField);

CMISSError CMISSEquationsSetSourceDestroyNum(const int RegionUserNumber,
		const int EquationsSetUserNumber);

CMISSError CMISSEquationsSetSourceDestroy(CMISSEquationsSetType *EquationsSet);

CMISSError CMISSEquationsSetSpecificationGetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *EquationsSetClass,
		int *EquationsSetType,
		int *EquationsSetSubtype);

CMISSError CMISSEquationsSetSpecificationGet(const CMISSEquationsSetType EquationsSet,
		int *EquationsSetClass,
		int *EquationsSetType,
		int *EquationsSetSubtype);

CMISSError CMISSEquationsSetSpecificationSetNum(const int RegionUserNumber,
		const int EquationsSetUserNumber,
		const int EquationsSetClass,
		const int EquationsSetType,
		const int EquationsSetSubtype);

CMISSError CMISSEquationsSetSpecificationSet(CMISSEquationsSetType *EquationsSet,
		const int EquationsSetClass,
		const int EquationsSetType,
		const int EquationsSetSubtype);

/*
 *==================================================================================================================================
 *
 * FIELD_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSFieldComponentInterpolationGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int ComponentNumber,
		int *InterpolationType);

CMISSError CMISSFieldComponentInterpolationGet(const CMISSFieldType Field,
		const int VariableType,
		const int ComponentNumber,
		int *InterpolationType);

CMISSError CMISSFieldComponentInterpolationSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int ComponentNumber,
		const int InterpolationType);

CMISSError CMISSFieldComponentInterpolationSet(const CMISSFieldType Field,
		const int VariableType,
		const int ComponentNumber,
		const int InterpolationType);

CMISSError CMISSFieldComponentLabelGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int ComponentNumber,
		const int LabelSize,
		char *Label);

CMISSError CMISSFieldComponentLabelGet(const CMISSFieldType Field,
		const int VariableType,
		const int ComponentNumber,
		const int LabelSize,
		char *Label);

CMISSError CMISSFieldComponentLabelSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int ComponentNumber,
		const int LabelSize,
		const char *Label);

CMISSError CMISSFieldComponentLabelSet(const CMISSFieldType Field,
		const int VariableType,
		const int ComponentNumber,
		const int LabelSize,
		const char *Label);

CMISSError CMISSFieldComponentMeshComponentGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
	    const int VariableType,
	    const int ComponentNumber,
	    int *MeshComponent);

CMISSError CMISSFieldComponentMeshComponentGet(const CMISSFieldType Field,
		const int VariableType,
		const int ComponentNumber,
		int *MeshComponent);

CMISSError CMISSFieldComponentMeshComponentSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int ComponentNumber,
		const int MeshComponent);

CMISSError CMISSFieldComponentMeshComponentSet(const CMISSFieldType Field,
		const int VariableType,
		const int ComponentNumber,
		const int MeshComponent);

CMISSError CMISSFieldComponentValuesInitialiseIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldComponentValuesInitialiseIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldComponentValuesInitialiseSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldComponentValuesInitialiseSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldComponentValuesInitialiseDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldComponentValuesInitialiseDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldComponentValuesInitialiseLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldComponentValuesInitialiseL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldDataTypeGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		int *DataType);

CMISSError CMISSFieldDataTypeGet(const CMISSFieldType Field,
		const int VariableType,
		int *DataType);

CMISSError CMISSFieldDataTypeSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int DataType);

CMISSError CMISSFieldDataTypeSet(const CMISSFieldType Field,
		const int VariableType,
		const int DataType);

CMISSError CMISSFieldDOFOrderTypeGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		int *DOFOrderType);

CMISSError CMISSFieldDOFOrderTypeGet(const CMISSFieldType Field,
		const int VariableType,
		int *DOFOrderType);

CMISSError CMISSFieldDOFOrderTypeSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int DOFOrderType);

CMISSError CMISSFieldDOFOrderTypeSet(const CMISSFieldType Field,
		const int VariableType,
		const int DOFOrderType);

CMISSError CMISSFieldCreateFinishNum(const int RegionUserNumber,
		const int FieldUserNumber);

CMISSError CMISSFieldCreateFinish(CMISSFieldType *Field);

CMISSError CMISSFieldCreateStartNum(const int RegionUserNumber,
		const int FieldUserNumber);

CMISSError CMISSFieldCreateStart(const int FieldUserNumber,
		const CMISSRegionType Region,
		const CMISSFieldType Field);

CMISSError CMISSFieldDependentTypeGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		int *DependentType);

CMISSError CMISSFieldDependentTypeGet(const CMISSFieldType Field,
		int *DependentType);

CMISSError CMISSFieldDependentTypeSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int DependentType);

CMISSError CMISSFieldDependentTypeSet(const CMISSFieldType Field,
		const int DependentType);

CMISSError CMISSFieldDestroyNum(const int RegionUserNumber,
		const int FieldUserNumber);

CMISSError CMISSFieldDestroy(const CMISSFieldType Field);

CMISSError CMISSFieldDimensionGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		int *Dimension);

CMISSError CMISSFieldDimensionGet(const CMISSFieldType Field,
		const int VariableType,
		int *Dimension);

CMISSError CMISSFieldDimensionSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int Dimension);

CMISSError CMISSFieldDimensionSet(const CMISSFieldType Field,
		const int VariableType,
		const int Dimension);

CMISSError CMISSFieldGeometricFieldGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		int *GeometricFieldUserNumber);

CMISSError CMISSFieldGeometricFieldGet(const CMISSFieldType Field,
		CMISSFieldType *GeometricField);

CMISSError CMISSFieldGeometricFieldSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int GeometricFieldUserNumber);

CMISSError CMISSFieldGeometricFieldSet(const CMISSFieldType Field,
		const CMISSFieldType GeometricField);

CMISSError CMISSFieldLabelGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int LabelSize,
		char *Label);

CMISSError CMISSFieldLabelGet(const CMISSFieldType Field,
		const int LabelSize,
		char *Label);

CMISSError CMISSFieldLabelSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int LabelSize,
		const char *Label);

CMISSError CMISSFieldLabelSet(const CMISSFieldType Field,
		const int LabelSize,
		const char *Label);

CMISSError CMISSFieldMeshDecompositionGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		int *DecompositionUserNumber);

CMISSError CMISSFieldMeshDecompositionGet(const CMISSFieldType Field,
		CMISSDecompositionType *MeshDecomposition);

CMISSError CMISSFieldMeshDecompositionSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber);

CMISSError CMISSFieldMeshDecompositionSet(const CMISSFieldType Field,
		const CMISSDecompositionType MeshDecomposition);

CMISSError CMISSFieldNumberOfComponentsGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		int *NumberOfComponents);

CMISSError CMISSFieldNumberOfComponentsGet(const CMISSFieldType Field,
		const int VariableType,
		int *NumberOfComponents);

CMISSError CMISSFieldNumberOfComponentsSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int NumberOfComponents);

CMISSError CMISSFieldNumberOfComponentsSet(const CMISSFieldType Field,
		const int VariableType,
		const int NumberOfComponents);

CMISSError CMISSFieldNumberOfVariablesGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		int *NumberOfVariables);

CMISSError CMISSFieldNumberOfVariablesGet(const CMISSFieldType Field,
		int *NumberOfVariables);

CMISSError CMISSFieldNumberOfVariablesSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int NumberOfVariables);

CMISSError CMISSFieldNumberOfVariablesSet(const CMISSFieldType Field,
		const int NumberOfVariables);

CMISSError CMISSFieldParameterSetAddConstantIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddConstantIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddConstantSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetAddConstantSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetAddConstantDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetAddConstantDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetAddConstantLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddConstantL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddElementIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddElementIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddElementSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetAddElementSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetAddElementDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetAddElementDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetAddElementLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddElementL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddNodeIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddNodeIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddNodeSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetAddNodeSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetAddNodeDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetAddNodeDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetAddNodeLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetAddNodeL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetCreateNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetCreate(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetDestroyNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetDestroy(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetDataGetIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataGetIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataGetSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		float *Parameters);

CMISSError CMISSFieldParameterSetDataGetSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		float *Parameters);

CMISSError CMISSFieldParameterSetDataGetDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		double *Parameters);

CMISSError CMISSFieldParameterSetDataGetDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int ParametersSize[1],
		double *Parameters);

CMISSError CMISSFieldParameterSetDataGetLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataGetL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		float *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		float *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		double *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		double *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetDataRestoreL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		int *ParametersSize[1],
		int *Parameters);

CMISSError CMISSFieldParameterSetGetConstantIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetConstantIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetConstantSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		float *Value);

CMISSError CMISSFieldParameterSetGetConstantSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		float *Value);

CMISSError CMISSFieldParameterSetGetConstantDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		double *Value);

CMISSError CMISSFieldParameterSetGetConstantDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		double *Value);

CMISSError CMISSFieldParameterSetGetConstantLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetConstantL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetElementIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementType,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetElementIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetElementSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		float *Value);

CMISSError CMISSFieldParameterSetGetElementSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		float *Value);

CMISSError CMISSFieldParameterSetGetElementDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		double *Value);

CMISSError CMISSFieldParameterSetGetElementDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		double *Value);

CMISSError CMISSFieldParameterSetGetElementLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		int Value);

CMISSError CMISSFieldParameterSetGetElementL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetNodeIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetNodeIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetNodeSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		float *Value);

CMISSError CMISSFieldParameterSetGetNodeSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		float *Value);

CMISSError CMISSFieldParameterSetGetNodeDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		double *Value);

CMISSError CMISSFieldParameterSetGetNodeDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		double *Value);

CMISSError CMISSFieldParameterSetGetNodeLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetGetNodeL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		int *Value);

CMISSError CMISSFieldParameterSetUpdateConstantIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateConstantIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateConstantSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetUpdateConstantSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetUpdateConstantDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetUpdateConstantDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetUpdateConstantLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateConstantL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateElementIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateElementIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateElementSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetUpdateElementSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetUpdateElementDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetUpdateElementDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetUpdateElementLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateElementL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int UserElementNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateFinishNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetUpdateFinish(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetUpdateNodeIntgNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateNodeIntg(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateNodeSPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetUpdateNodeSP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const float Value);

CMISSError CMISSFieldParameterSetUpdateNodeDPNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetUpdateNodeDP(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const double Value);

CMISSError CMISSFieldParameterSetUpdateNodeLNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateNodeL(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType,
		const int DerivativeNumber,
		const int UserNodeNumber,
		const int ComponentNumber,
		const int Value);

CMISSError CMISSFieldParameterSetUpdateStartNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParameterSetUpdateStart(const CMISSFieldType Field,
		const int VariableType,
		const int FieldSetType);

CMISSError CMISSFieldParametersToFieldParametersComponentCopyNum(const int FromRegionUserNumber,
		const int FromFieldUserNumber,
		const int FromVariableType,
		const int FromParameterSetType,
		const int FromComponentNumber,
		const int ToRegionUserNumber,
		const int ToFieldUserNumber,
		const int ToVariableType,
		const int ToParameterSetType,
		const int ToComponentNumber);

CMISSError CMISSFieldParametersToFieldParametersComponentCopy(const CMISSFieldType FromField,
		const int FromVariableType,
		const int FromParameterSetType,
		const int FromComponentNumber,
		const CMISSFieldType ToField,
		const int ToVariableType,
		const int ToParameterSetType,
		const int ToComponentNumber);

CMISSError CMISSFieldScalingTypeGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		int *ScalingType);

CMISSError CMISSFieldScalingTypeGet(const CMISSFieldType Field,
		int *ScalingType);

CMISSError CMISSFieldScalingTypeSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int ScalingType);

CMISSError CMISSFieldScalingTypeSet(const CMISSFieldType Field,
		const int ScalingType);

CMISSError CMISSFieldTypeGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		int *FieldType);

CMISSError CMISSFieldTypeGet(const CMISSFieldType Field,
		int *FieldType);

CMISSError CMISSFieldTypeSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int FieldType);

CMISSError CMISSFieldTypeSet(const CMISSFieldType Field,
		const int FieldType);

CMISSError CMISSFieldVariableLabelGetCNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int LabelSize,
		char *Label);

CMISSError CMISSFieldVariableLabelGetC(const CMISSFieldType Field,
		const int VariableType,
		const int LabelSize,
		char *Label);

CMISSError CMISSFieldVariableLabelSetCNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableType,
		const int LabelSize,
		const char *Label);

CMISSError CMISSFieldVariableLabelSetC(const CMISSFieldType Field,
		const int VariableType,
		const int LabelSize,
		const char *Label);

CMISSError CMISSFieldVariableTypesGetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		int *VariableTypesSize[1],
		int *VariableTypes);

CMISSError CMISSFieldVariableTypesGet(const CMISSFieldType Field,
		int *VariableTypesSize[1],
		int *VariableTypes);

CMISSError CMISSFieldVariableTypesSetNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int VariableTypesSize[1],
		const int *VariableTypes);

CMISSError CMISSFieldVariableTypesSet(const CMISSFieldType Field,
		const int VariableTypesSize[1],
		const int *VariableTypes);

/*
 *==================================================================================================================================
 *
 * FIELD_IO_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSFieldIOElementsExportC(CMISSFieldsType *Fields,
		const int FileNameSize,
		const char *FileName,
		const int MethodSize,
		const char *Method);

CMISSError CMISSFieldIONodesExportC(CMISSFieldsType *Fields,
		const int FileNameSize,
		const char *FileName,
		const int MethodSize,
		const char *Method);\

/*
 *==================================================================================================================================
 *
 * GENERATED_MESH_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSGeneratedMeshBasisGetNum(const int GeneratedMeshUserNumber,
		int *BasisUserNumber);

CMISSError CMISSGeneratedMeshBasisGet(const CMISSGeneratedMeshType GeneratedMesh,
		CMISSBasisType *Basis);

CMISSError CMISSGeneratedMeshBasisSetNum(const int GeneratedMeshUserNumber,
		const int BasisUserNumber);

CMISSError CMISSGeneratedMeshBasisSet(const CMISSGeneratedMeshType GeneratedMesh,
		const CMISSBasisType Basis);

CMISSError CMISSGeneratedMeshCreateFinishNum(const int GeneratedMeshUserNumber,
		const int BasisUserNumber);

CMISSError CMISSGeneratedMeshCreateFinish(const CMISSGeneratedMeshType GeneratedMesh,
		const int MeshUserNumber,
		CMISSMeshType *Mesh);

CMISSError CMISSGeneratedMeshCreateStartNum(const int GeneratedMeshUserNumber,
		const int RegionUserNumber);

CMISSError CMISSGeneratedMeshCreateStart(const int GeneratedMeshUserNumber,
		CMISSRegionType *Region,
		CMISSGeneratedMeshType *GeneratedMesh);

CMISSError CMISSGeneratedMeshDestroyNum(const int GeneratedMeshUserNumber);

CMISSError CMISSGeneratedMeshDestroy(CMISSGeneratedMeshType *GeneratedMesh);

CMISSError CMISSGeneratedMeshExtentGetNum(const int GeneratedMeshUserNumber,
		int *ExtentSize,
		int *Extent);

CMISSError CMISSGeneratedMeshExtentGet(const CMISSGeneratedMeshType GeneratedMesh,
		int *ExtentSize,
		int *Extent);

CMISSError CMISSGeneratedMeshExtentSetNum(const int GeneratedMeshUserNumber,
		const int ExtentSize[1],
		const int *Extent);

CMISSError CMISSGeneratedMeshExtentSet(const CMISSGeneratedMeshType GeneratedMesh,
		const int ExtentSize[1],
		const int *Extent);

CMISSError CMISSGeneratedMeshNumberOfElementsGetNum(const int GeneratedMeshUserNumber,
		int *ExtentSize,
		int *Extent);

CMISSError CMISSGeneratedMeshNumberOfElementsGet(const CMISSGeneratedMeshType GeneratedMesh,
		int *ExtentSize,
		int *Extent);

CMISSError CMISSGeneratedMeshNumberOfElementsSetNum(const int GeneratedMeshUserNumber,
		const int ExtentSize[1],
		const int *Extent);

CMISSError CMISSGeneratedMeshNumberOfElementsSet(const CMISSGeneratedMeshType GeneratedMesh,
		const int ExtentSize[1],
		const int *Extent);

CMISSError CMISSGeneratedMeshOriginGetNum(const int GeneratedMeshUserNumber,
		int *OriginSize,
		int *Origin);

CMISSError CMISSGeneratedMeshOriginGet(const CMISSGeneratedMeshType GeneratedMesh,
		int *OriginSize,
		int *Origin);

CMISSError CMISSGeneratedMeshOriginSetNum(const int GeneratedMeshUserNumber,
		const int OriginSize[1],
		const int *Origin);

CMISSError CMISSGeneratedMeshOriginSet(const CMISSGeneratedMeshType GeneratedMesh,
		const int OriginSize[1],
		const int *Origin);

CMISSError CMISSGeneratedMeshTypeGetNum(const int GeneratedMeshUserNumber,
		int *GeneratedMeshType);

CMISSError CMISSGeneratedMeshTypeGet(const CMISSGeneratedMeshType GeneratedMesh,
		int *GeneratedMeshType);

CMISSError CMISSGeneratedMeshTypeSetNum(const int GeneratedMeshUserNumber,
		const int GeneratedMeshType);

CMISSError CMISSGeneratedMeshTypeSet(const CMISSGeneratedMeshType GeneratedMesh,
		const int GeneratedMeshType);

CMISSError CMISSGeneratedMeshGeometricParametersCalculateNum(const int RegionUserNumber,
		const int FieldUserNumber,
		const int GeneratedMeshUserNumber);

CMISSError CMISSGeneratedMeshGeometricParametersCalculate(CMISSFieldType *Field,
		const CMISSGeneratedMeshType GeneratedMesh);

/*
 *==================================================================================================================================
 *
 * MESH_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSDecompositionCreateFinishNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber);

CMISSError CMISSDecompositionCreateFinish(const CMISSDecompositionType Decomposition);

CMISSError CMISSDecompositionCreateStartNum(const int DecompositionUserNumber,
		const int RegionUserNumber,
		const int MeshUserNumber);

CMISSError CMISSDecompositionCreateStart(const int DecompositionUserNumber,
		const CMISSMeshType Mesh,
		CMISSDecompositionType *Decomposition);

CMISSError CMISSDecompositionDestroyNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber);

CMISSError CMISSDecompositionDestroy(const CMISSDecompositionType Decomposition);

CMISSError CMISSDecompositionElementDomainCalculateNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber);

CMISSError CMISSDecompositionElementDomainCalculate(const CMISSDecompositionType Decomposition);

CMISSError CMISSDecompositionElementDomainGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		const int ElementUserNumber,
		int *Domain);

CMISSError CMISSDecompositionElementDomainGet(const CMISSDecompositionType Decomposition,
		const int ElementUserNumber,
		int *Domain);

CMISSError CMISSDecompositionElementDomainSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		const int ElementUserNumber,
		const int Domain);

CMISSError CMISSDecompositionElementDomainSet(const CMISSDecompositionType Decomposition,
		const int ElementUserNumber,
		const int Domain);

CMISSError CMISSDecompositionMeshComponentGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		int *MeshComponentNumber);

CMISSError CMISSDecompositionMeshComponentGet(const CMISSDecompositionType Decomposition,
		int *MeshComponentNumber);

CMISSError CMISSDecompositionMeshComponentSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		const int MeshComponentNumber);

CMISSError CMISSDecompositionMeshComponentSet(const CMISSDecompositionType Decomposition,
		const int MeshComponentNumber);

CMISSError CMISSDecompositionNumberOfDomainsGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		int *NumberOfDomains);

CMISSError CMISSDecompositionNumberOfDomainsGet(const CMISSDecompositionType Decomposition,
		int *NumberOfDomains);

CMISSError CMISSDecompositionNumberOfDomainsSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		const int NumberOfDomains);

CMISSError CMISSDecompositionNumberOfDomainsSet(const CMISSDecompositionType Decomposition,
		const int NumberOfDomains);

CMISSError CMISSDecompositionTypeGetCNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		int *DecompositionType);

CMISSError CMISSDecompositionTypeGet(const CMISSDecompositionType Decomposition,
		int *DecompositionType);

CMISSError CMISSDecompositionTypeSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		const int DecompositionType);

CMISSError CMISSDecompositionTypeSet(const CMISSDecompositionType Decomposition,
		const int DecompositionType);

CMISSError CMISSDecompositionNodeDomainGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int DecompositionUserNumber,
		const int NodeUserNumber,
		const int MeshComponentNumber,
		int *Domain);

CMISSError CMISSDecompositionNodeDomainGet(const CMISSDecompositionType Decomposition,
		const int NodeUserNumber,
		const int MeshUserNumber,
		int *Domain);

CMISSError CMISSMeshCreateFinishNum(const int RegionUserNumber,
		const int MeshUserNumber);

CMISSError CMISSMeshCreateFinish(const CMISSMeshType Mesh);

CMISSError CMISSMeshCreateStartNum(const int MeshUserNumber,
		const int RegionUserNumber,
		const int NumberOfDimensions);

CMISSError CMISSMeshCreateStart(const int MeshUserNumber,
		const CMISSRegionType Region,
		const int NumberOfDimensions,
		CMISSMeshType *Mesh);

CMISSError CMISSMeshDestroyNum(const int RegionUserNumber,
		const int MeshUserNumber);

CMISSError CMISSMeshDestroy(const CMISSMeshType Mesh);

CMISSError CMISSMeshNumberOfComponentsGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		int *NumberOfComponents);

CMISSError CMISSMeshNumberOfComponentsGet(const CMISSMeshType Mesh,
		int *NumberOfComponents);

CMISSError CMISSMeshNumberOfComponentsSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int NumberOfComponents);

CMISSError CMISSMeshNumberOfComponentsSet(const CMISSMeshType Mesh,
		const int NumberOfComponents);

CMISSError CMISSMeshNumberOfElementsGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		int *NumberOfElements);

CMISSError CMISSMeshNumberOfElementsGet(const CMISSMeshType Mesh,
		int *NumberOfElements);

CMISSError CMISSMeshNumberOfElementsSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int NumberOfElements);

CMISSError CMISSMeshNumberOfElementsSet(const CMISSMeshType Mesh,
		const int NumberOfElements);

CMISSError CMISSMeshElementsCreateFinishNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber);

CMISSError CMISSMeshElementsCreateFinish(const CMISSMeshElementsType MeshElements);

CMISSError CMISSMeshElementsCreateStartNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int BasisUserNumber);

CMISSError CMISSMeshElementsCreateStart(CMISSMeshType *Mesh,
		const int MeshComponentNumber,
		const CMISSBasisType Basis,
		CMISSMeshElementsType *MeshElements);

CMISSError CMISSMeshElementsBasisGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int GlobalElementNumber,
		int *BasisUserNumber);

CMISSError CMISSMeshElementsBasisGet(const CMISSMeshElementsType MeshElements,
		const int GlobalElementNumber,
		CMISSBasisType *Basis);

CMISSError CMISSMeshElementsBasisSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int GlobalElementNumber,
		const int BasisUserNumber);

CMISSError CMISSMeshElementsBasisSet(const CMISSMeshElementsType MeshElements,
		const int GlobalElementNumber,
		const CMISSBasisType Basis);

CMISSError CMISSMeshElementsNodesGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int GlobalElementNumber,
		int *ElementUserNodesSize,
		int *ElementUserNodes);

CMISSError CMISSMeshElementsNodesGet(const CMISSMeshElementsType MeshElements,
		const int GlobalElementNumber,
		int *ElementUserNodesSize,
		int *ElementUserNodes);

CMISSError CMISSMeshElementsNodesSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int GlobalElementNumber,
		const int ElementUserNodesSize[1],
		const int *ElementUserNodes);

CMISSError CMISSMeshElementsNodesSet(const CMISSMeshElementsType MeshElements,
		const int GlobalElementNumber,
		const int ElementUserNodesSize[1],
		const int *ElementUserNodes);

CMISSError CMISSMeshElementsUserNumberGetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int ElementGlobalNumber,
		int *ElementUserNumber);

CMISSError CMISSMeshElementsUserNumberGet(const CMISSMeshElementsType MeshElements,
		const int ElementGlobalNumber,
		int *ElementUserNumber);

CMISSError CMISSMeshElementsUserNumberSetNum(const int RegionUserNumber,
		const int MeshUserNumber,
		const int MeshComponentNumber,
		const int ElementGlobalNumber,
		const int ElementUserNumber);

CMISSError CMISSMeshElementsUserNumberSet(const CMISSMeshElementsType MeshElements,
		const int ElementGlobalNumber,
		const int ElementUserNumber);

/*
 *==================================================================================================================================
 *
 * NODE_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSNodesCreateFinishNum(const int RegionUserNumber);

CMISSError CMISSNodesCreateFinish(const CMISSNodesType Nodes);

CMISSError CMISSNodesCreateStartNum(const int RegionUserNumber,
		const int NumberOfNodes);

CMISSError CMISSNodesCreateStart(const CMISSRegionType Region,
		const int NumberOfNodes,
		const CMISSNodesType Nodes);

CMISSError CMISSNodesDestroyNum(const int RegionUserNumber);

CMISSError CMISSNodesDestroy(const CMISSNodesType Nodes);

CMISSError CMISSNodesLabelGetCNum(const int RegionUserNumber,
		const int NodeGlobalNumber,
		const int LabelSize,
		char *Label);

CMISSError CMISSNodesLabelGetC(const CMISSNodesType Nodes,
		const int NodeGlobalNumber,
		const int LabelSize,
		char *Label);

CMISSError CMISSNodesLabelSetCNum(const int RegionUserNumber,
		const int NodeGlobalNumber,
		const int LabelSize,
		const char *Label);

CMISSError CMISSNodesLabelSetC(const CMISSNodesType Nodes,
		const int NodeGlobalNumber,
		const int LabelSize,
		const char *Label);

CMISSError CMISSNodesUserNumberGetNum(const int RegionUserNumber,
		const int NodeGlobalNumber,
		int *NodeUserNumber);

CMISSError CMISSNodesUserNumberGet(const CMISSNodesType Nodes,
		const int NodeGlobalNumber,
		int *NodeUserNumber);

CMISSError CMISSNodesUserNumberSetNum(const int RegionUserNumber,
		const int NodeGlobalNumber,
		const int NodeUserNumber);

CMISSError CMISSNodesUserNumberSet(const CMISSNodesType Nodes,
		const int NodeGlobalNumber,
		const int NodeUserNumber);

/*
 *==================================================================================================================================
 *
 * PROBLEM_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSProblemCreateFinishNum(const int ProblemUserNumber);

CMISSError CMISSProblemCreateFinish(const CMISSProblemType Problem);

CMISSError CMISSProblemCreateStartNum(const int ProblemUserNumber);

CMISSError CMISSProblemCreateStart(const int ProblemUserNumber,
		const CMISSProblemType Problem);

CMISSError CMISSProblemDestroyNum(const int ProblemUserNumber);

CMISSError CMISSProblemDestroy(const CMISSProblemType Problem);

CMISSError CMISSProblemControlLoopCreateFinishNum(const int ProblemUserNumber);

CMISSError CMISSProblemControlLoopCreateFinish(const CMISSProblemType Problem);

CMISSError CMISSProblemControlLoopCreateStartNum(const int ProblemUserNumber);

CMISSError CMISSProblemControlLoopCreateStart(const CMISSProblemType Problem);

CMISSError CMISSProblemControlLoopDestroyNum(const int ProblemUserNumber);

CMISSError CMISSProblemControlLoopDestroy(const CMISSProblemType Problem);

CMISSError CMISSProblemControlLoopGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		CMISSControlLoopType *ControlLoop);

CMISSError CMISSProblemControlLoopGet(const CMISSProblemType Problem,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		CMISSControlLoopType *ControlLoop);

CMISSError CMISSProblemSolveNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolve(const CMISSProblemType Problem);

CMISSError CMISSProblemSolverGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		CMISSSolverType *Solver);

CMISSError CMISSProblemSolverGet(const CMISSProblemType Problem,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		CMISSSolverType *Solver);

CMISSError CMISSProblemSolverEquationsCreateFinishNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolverEquationsCreateFinish(const CMISSProblemType Problem);

CMISSError CMISSProblemSolverEquationsCreateStartNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolverEquationsCreateStart(const CMISSProblemType Problem);

CMISSError CMISSProblemSolverEquationsDestroyNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolverEquationsDestroy(const CMISSProblemType Problem);

CMISSError CMISSProblemSolverEquationsGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		CMISSSolverEquationsType *SolverEquations);

CMISSError CMISSProblemSolverEquationsGet(const CMISSProblemType Problem,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		CMISSSolverEquationsType *SolverEquations);

CMISSError CMISSProblemSolversCreateFinishNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolversCreateFinish(const CMISSProblemType Problem);

CMISSError CMISSProblemSolversCreateStartNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolversCreateStart(const CMISSProblemType Problem);

CMISSError CMISSProblemSolversDestroyNum(const int ProblemUserNumber);

CMISSError CMISSProblemSolversDestroy(const CMISSProblemType Problem);

CMISSError CMISSProblemSpecificationGetNum(const int ProblemUserNumber,
		int *ProblemClass,
		int *ProblemType,
		int *ProblemSubtype);

CMISSError CMISSProblemSpecificationGet(const CMISSProblemType Problem,
		int *ProblemClass,
		int *ProblemType,
		int *ProblemSubtype);

CMISSError CMISSProblemSpecificationSetNum(const int ProblemUserNumber,
		const int ProblemClass,
		const int ProblemType,
		const int ProblemSubtype);

CMISSError CMISSProblemSpecificationSet(const CMISSProblemType Problem,
		const int ProblemClass,
		const int ProblemType,
		const int ProblemSubtype);

/*
 *=================================================================================================================================
 *
 * REGION
 *
 *==================================================================================================================================
 */

CMISSError CMISSRegionCoordinateSystemGetNum(const int RegionUserNumber,
		int *CoordinateSystemUserNumber);

CMISSError CMISSRegionCoordinateSystemGet(const CMISSRegionType Region,
		CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSRegionCoordinateSystemSetNum(const int RegionUserNumber,
		const int CoordinateSystemUserNumber);

CMISSError CMISSRegionCoordinateSystemSet(const CMISSRegionType Region,
		const CMISSCoordinateSystemType CoordinateSystem);

CMISSError CMISSRegionCreateFinishNum(const int RegionUserNumber);

CMISSError CMISSRegionCreateFinish(const CMISSRegionType Region);

CMISSError CMISSRegionCreateStartNum(const int RegionUserNumber,
				     const int ParentRegionUserNumber);

CMISSError CMISSRegionCreateStart(const int RegionUserNumber,
				  const CMISSRegionType ParentRegion,
				  CMISSRegionType *Region);

CMISSError CMISSRegionDestroyNum(const int RegionUserNumber);

CMISSError CMISSRegionDestroy(CMISSRegionType *Region);

CMISSError CMISSRegionLabelGetNum(const int RegionUserNumber,
				  const int LabelSize,
				  char *Label);

CMISSError CMISSRegionLabelGet(const CMISSRegionType Region,
			       const int LabelSize,
			       char *Label);

CMISSError CMISSRegionLabelSetNum(const int RegionUserNumber,
				  const int LabelSize,
				  const char *Label);

CMISSError CMISSRegionLabelSet(const CMISSRegionType Region,
			       const int LabelSize,
			       const char *Label);

/*
 *==================================================================================================================================
 *
 * SOLVER_ROUTINES
 *
 *==================================================================================================================================
 */

CMISSError CMISSSolverDAEEulerSolverTypeGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *DAEEulerSolverType);

CMISSError CMISSSolverDAEEulerSolverTypeGet(const CMISSSolverType Solver,
		int *DAEEulerSolverType);

CMISSError CMISSSolverDAEEulerSolverTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int DAEEulerSolverType);

CMISSError CMISSSolverDAEEulerSolverTypeSet(const CMISSSolverType Solver,
		const int DAEEulerSolverType);

CMISSError CMISSSolverDAESolverTypeGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *DAESolverType);

CMISSError CMISSSolverDAESolverTypeGet(const CMISSSolverType Solver,
		int *DAESolverType);

CMISSError CMISSSolverDAESolverTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int DAESolverType);

CMISSError CMISSSolverDAESolverTypeSet(const CMISSSolverType Solver,
		const int DAESolverType);

CMISSError CMISSSolverDAETimesSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int StartTime,
		const int EndTime,
		const int IntialStep);

CMISSError CMISSSolverDAETimesSet(const CMISSSolverType Solver,
		const int StartTime,
		const int EndTime,
		const int InitialStep);

CMISSError CMISSSolverDynamicDegreeGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *Degree);

CMISSError CMISSSolverDynamicDegreeGet(const CMISSSolverType Solver,
		int *Degree);

CMISSError CMISSSolverDynamicDegreeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int Degree);

CMISSError CMISSSolverDynamicDegreeSet(const CMISSSolverType Solver,
		const int Degree);

CMISSError CMISSSolverDynamicLinearityTypeGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *LinearityType);

CMISSError CMISSSolverDynamicLinearityTypeGet(const CMISSSolverType Solver,
		int *LinearityType);

CMISSError CMISSSolverDynamicNonlinearSolverGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *NonlinearSolverIndex);

CMISSError CMISSSolverDynamicNonlinearSolverGet(const CMISSSolverType Solver,
		CMISSSolverType *NonLinearSolver);

CMISSError CMISSSolverDynamicLinearSolverGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *LinearSolverIndex);

CMISSError CMISSSolverDynamicLinearSolverGet(const CMISSSolverType Solver,
		CMISSSolverType *LinearSolver);

CMISSError CMISSSolverDynamicSchemeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int Scheme);

CMISSError CMISSSolverDynamicSchemeSet(const CMISSSolverType Solver,
		const int Scheme);

CMISSError CMISSSolverDynamicThetaSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int ThetasSize[1],
		const double *Thetas);

CMISSError CMISSSolverDynamicThetaSet(const CMISSSolverType Solver,
		const int ThetasSize[1],
		const double *Thetas);

CMISSError CMISSSolverDynamicTimesSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double CurrentTime,
		const double TimeIncrement);

CMISSError CMISSSolverDynamicTimesSet(const CMISSSolverType Solver,
		const double CurrentTime,
		const double TimeIncrement);

CMISSError CMISSSolverLibraryTypeGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *LibraryType);

CMISSError CMISSSolverLibraryTypeGet(const CMISSSolverType Solver,
		int *LibraryType);

CMISSError CMISSSolverLibraryTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int LibraryType);

CMISSError CMISSSolverLibraryTypeSet(const CMISSSolverType Solver,
		const int LibraryType);

CMISSError CMISSSolverLinearDirectTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int DirectSolverType);

CMISSError CMISSSolverLinearDirectTypeSet(const CMISSSolverType Solver,
		const int DirectSolverType);

CMISSError CMISSSolverLinearIterativeAbsoluteToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double AbsoluteTolerance);

CMISSError CMISSSolverLinearIterativeAbsoluteToleranceSet(const CMISSSolverType Solver,
		const double AbsoluteTolerance);

CMISSError CMISSSolverLinearIterativeDivergenceToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double DivergenceTolerance);

CMISSError CMISSSolverLinearIterativeDivergenceToleranceSet(const CMISSSolverType Solver,
		const double DivergenceTolerance);

CMISSError CMISSSolverLinearIterativeGMRESRestartSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int GMRESRestart);

CMISSError CMISSSolverLinearIterativeGMRESRestartSet(const CMISSSolverType Solver,
		const int GMRESRestart);

CMISSError CMISSSolverLinearIterativeMaximumIterationsSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int MaximumIterations);

CMISSError CMISSSolverLinearIterativeMaximumIterationsSet(const CMISSSolverType Solver,
		const int MaximumIterations);

CMISSError CMISSSolverLinearIterativePreconditionerTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int PreconditionerType);

CMISSError CMISSSolverLinearIterativePreconditionerTypeSet(const CMISSSolverType Solver,
		const int PreconditionerType);

CMISSError CMISSSolverLinearIterativeRelativeToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double RelativeTolerance);

CMISSError CMISSSolverLinearIterativeRelativeToleranceSet(const CMISSSolverType Solver,
		const double RelativeTolerance);

CMISSError CMISSSolverLinearIterativeTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int IterativeSolverType);

CMISSError CMISSSolverLinearIterativeTypeSet(const CMISSSolverType Solver,
		const int IterativeSolverType);

CMISSError CMISSSolverLinearTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int LinearSolverType);

CMISSError CMISSSolverLinearTypeSet(const CMISSSolverType Solver,
		const int LinearSolverType);

CMISSError CMISSSolverNewtonAbsoluteToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double AbsoluteTolerance);

CMISSError CMISSSolverNewtonAbsoluteToleranceSet(const CMISSSolverType Solver,
		const double AbsoluteTolerance);

CMISSError CMISSSolverNewtonJacobianCalculationTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int JacobianCalculationType);

CMISSError CMISSSolverNewtonJacobianCalculationTypeSet(const CMISSSolverType Solver,
		const int JacobianCalculationType);

CMISSError CMISSSolverNewtonLinearSolverGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		int *LinearSolverIndex);

CMISSError CMISSSolverNewtonLinearSolverGet(const CMISSSolverType Solver,
		CMISSSolverType *LinearSolver);

CMISSError CMISSSolverNewtonLineSearchAlphaSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double Alpha);

CMISSError CMISSSolverNewtonLineSearchAlphaSet(const CMISSSolverType Solver,
		const double Alpha);

CMISSError CMISSSolverNewtonLineSearchMaxStepSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double MaxStep);

CMISSError CMISSSolverNewtonLineSearchMaxStepSet(const CMISSSolverType Solver,
		const double MaxStep);

CMISSError CMISSSolverNewtonLineSearchStepTolSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double StepTol);

CMISSError CMISSSolverNewtonLineSearchStepTolSet(const CMISSSolverType Solver,
		const double StepTol);

CMISSError CMISSSolverNewtonLineSearchTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int LineSearchType);

CMISSError CMISSSolverNewtonLineSearchTypeSet(const CMISSSolverType Solver,
		const int LineSearchType);

CMISSError CMISSSolverNewtonMaximumFunctionEvaluationsSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int MaximumFunctionEvaluations);

CMISSError CMISSSolverNewtonMaximumFunctionEvaluationsSet(const CMISSSolverType Solver,
		const int MaximumFunctionEvaluations);

CMISSError CMISSSolverNewtonMaximumIterationsSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int MaximumIterations);

CMISSError CMISSSolverNewtonMaximumIterationsSet(const CMISSSolverType Solver,
		const int MaximumIterations);

CMISSError CMISSSolverNewtonRelativeToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double RelativeTolerance);

CMISSError CMISSSolverNewtonRelativeToleranceSet(const CMISSSolverType Solver,
		const double RelativeTolerance);

CMISSError CMISSSolverNewtonSolutionToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double SolutionTolerance);

CMISSError CMISSSolverNewtonSolutionToleranceSet(const CMISSSolverType Solver,
		const double SolutionTolerance);

CMISSError CMISSSolverNewtonTrustRegionDelta0SetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const double Delta0);

CMISSError CMISSSolverNewtonTrustRegionDelta0Set(const CMISSSolverType Solver,
		const double Delta0);

CMISSError CMISSSolverNewtonTrustRegionToleranceSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int Tolerance);

CMISSError CMISSSolverNewtonTrustRegionToleranceSet(const CMISSSolverType Solver,
		const int Tolerance);

CMISSError CMISSSolverNewtonTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int NewtonSolveType);

CMISSError CMISSSolverNewtonTypeSet(const CMISSSolverType Solver,
		const int NewtonSolveType);

CMISSError CMISSSolverNonlinearTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int NonlinearSolveType);

CMISSError CMISSSolverNonlinearTypeSet(const CMISSSolverType Solver,
		const int NonlinearSolveType);

CMISSError CMISSSolverOutputTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int OutputType);

CMISSError CMISSSolverOutputTypeSet(const CMISSSolverType Solver,
		const int OutputType);

CMISSError CMISSSolverSolverEquationsGetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		CMISSSolverEquationsType *SolverEquations);

CMISSError CMISSSolverSolverEquationsGet(const CMISSSolverType Solver,
		CMISSSolverEquationsType *SolverEquations);

CMISSError CMISSSolverEquationsEquationsSetAddNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int RegionUserNumber,
		const int EquationsSetUserNumber,
		int *EquationsSetIndex);

CMISSError CMISSSolverEquationsEquationsSetAdd(const CMISSSolverEquationsType SolverEquations,
		const CMISSEquationsSetType EquationsSet,
		int *EquationsSetIndex);

CMISSError CMISSSolverEquationsSparsityTypeSetNum(const int ProblemUserNumber,
		const int ControlLoopIdentifiersSize[1],
		const int *ControlLoopIdentifiers,
		const int SolverIndex,
		const int SparsityType);

CMISSError CMISSSolverEquationsSparsityTypeSet(const CMISSSolverEquationsType SolverEquations,
		const int SparsityType);
