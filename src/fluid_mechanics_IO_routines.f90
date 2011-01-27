!> \file
!> $Id$
!> \author Sebastian Krittian
!> \brief This module handles some mesh/parameter input routines and cmgui output routines for fluid mechanics
!> routines and should be eventually replaces by field_IO_routines.f90 
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
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
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

!> Temporary IO routines for fluid mechanics

MODULE FLUID_MECHANICS_IO_ROUTINES

 USE BASE_ROUTINES
 USE EQUATIONS_SET_CONSTANTS
 USE FIELD_ROUTINES
 USE TYPES
 USE INPUT_OUTPUT 
 USE KINDS   

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !1=M, 2=V, 3=P !
  TYPE ARRAY_MESH
    INTEGER(INTG) Lx                                    ! External Structure: Number of Tessellation Coefficients
    INTEGER(INTG) Lt                                    ! External Structure: Number of Tessellation elements
    INTEGER(INTG) ID                                    ! Internal Structure: Basis number to call letter map
    INTEGER(INTG) Lb
    REAL(DP), pointer :: X(:,:)                         ! External Structure: Tessellation coefficients
    INTEGER(INTG), pointer :: T(:,:)                    ! External Structure: Tessellation map
    INTEGER(INTG), pointer :: NBT(:,:)                  ! Internal Structure: Neighboring element map
    INTEGER(INTG), pointer :: B(:,:)
  END TYPE ARRAY_MESH

  TYPE ARRAY_INT
    REAL(DP), pointer :: Y(:,:)                         ! Internal Structure: Basis @ given volume quadrature points
    REAL(DP), pointer :: Y_f(:,:,:)                     ! Internal Structure: Basis @ given facet quadrature points
    REAL(DP), pointer :: dY(:,:,:), TdY(:,:,:), dY_f(:,:,:,:)
  END TYPE ARRAY_INT

  TYPE ARRAY_BASE
    INTEGER(INTG) n                                     ! External Structure: Function Space dimension on T_e
    INTEGER(INTG) nl                                    ! Internal Structure: Breakdown for tensor input
    INTEGER(INTG) DISCONT                               ! External Structure: Continuity constraint
    INTEGER(INTG) DM                                    ! External Structure: Field dimension
    TYPE(ARRAY_INT) I                                   ! External Structure: Field integrations
    REAL(DP), pointer :: Q(:,:)                         ! External Structure: Basis coefficient tensor
    REAL(DP), pointer :: P(:,:)                         ! External Structure: Basis power tensor
    REAL(DP), pointer :: XI(:,:)                        ! External Structure: Basis xi-point coordinates
    INTEGER(INTG), pointer :: B_ID(:,:)                 ! External Structure: Basis ordering tensor
    CHARACTER*2 CL                                      ! External Structure: Basis call letter
  END TYPE ARRAY_BASE

  TYPE ARRAY_PROBLEM_BASE
    TYPE(ARRAY_BASE), pointer :: B(:)                   ! Internal Structure: Basis
    INTEGER(INTG) n_pts                                 ! External Structure: number of volume quadrature points
    INTEGER(INTG) n_pts_f                               ! External Structure: number of facet quadrature points
    INTEGER(INTG) n_ptsl                                ! Internal Structure: volume quadrature for tensor input
    INTEGER(INTG) n_pts_fl                              ! Internal Structure: facet quadrature for tensor input
    INTEGER(INTG) HEXA                                  ! External Structure: tensor input parameter
    INTEGER(INTG) FACES                                 ! External Structure: number of master element facets
    INTEGER(INTG) FNODES                                ! External Structure: number of nodes per facet (spatial map)
    INTEGER(INTG) n_B                                   ! External Structure: number of basis
    INTEGER(INTG) DM                                    ! External Structure: spatial dimension of master element
    INTEGER(INTG) SPL, VEL, PRS
    INTEGER(INTG) TRI_BASIS, QUAD_BASIS, TET_BASIS, HEX_BASIS
    REAL(DP) VL                                         ! External Structure: volume of master element
    REAL(DP) VL_f(6)                                    ! External Structure: facet volume of master element
    REAL(DP) nrm(3,6)                                   ! Internal Structure: facet normals of master element
    REAL(DP), pointer :: gpt(:,:)                       ! External Structure: volume quadrature points
    REAL(DP), pointer :: gw(:)                          ! External Structure: volume quadrature weights
    REAL(DP), pointer :: gpt_f(:,:,:)                   ! External Structure: facet quadrature points
    REAL(DP), pointer :: gw_f(:)                        ! External Structure: facet quadrature weights
  END TYPE ARRAY_PROBLEM_BASE

  TYPE EXPORT_CONTAINER
    REAL(DP), POINTER::N(:,:)
    INTEGER(INTG), POINTER::M(:,:),V(:,:),P(:,:)
    INTEGER(INTG):: D,F,ID_M,ID_V,ID_P,IT_M,IT_V,IT_P,IT_T
    INTEGER(INTG):: E_M,E_P,E_V,E_T,EN_M,EN_P,EN_V,EN_T,N_M,N_P,N_V,N_T
  END TYPE EXPORT_CONTAINER

  TYPE DARCY_PARAMETERS
    !chrm, 20.08.09: For passing data of Darcy problems
    INTEGER:: TESTCASE
    INTEGER:: BC_NUMBER_OF_WALL_NODES, NUMBER_OF_BCS
    REAL(DP):: PERM, VIS, PERM_OVER_VIS, P_SINK
    REAL(DP):: X1, X2, Y1, Y2, Z1, Z2
    REAL(DP):: LENGTH, GEOM_TOL
    REAL(DP):: max_node_spacing
    LOGICAL :: STAB, DEBUG, ANALYTIC
  END TYPE DARCY_PARAMETERS

  TYPE COUPLING_PARAMETERS
    INTEGER:: NUMBER_OF_COUPLINGS
    INTEGER, POINTER:: INTERFACE_ELEMENT_NUMBER(:)
    INTEGER, POINTER:: INTERFACE_ELEMENT_LOCAL_NODE(:)
    INTEGER:: MESH1_ID
    INTEGER, POINTER:: MESH1_ELEMENT_NUMBER(:)  
    REAL(DP), POINTER:: MESH1_ELEMENT_XI(:,:)
    INTEGER:: MESH2_ID
    INTEGER, POINTER:: MESH2_ELEMENT_NUMBER(:)  
    REAL(DP), POINTER:: MESH2_ELEMENT_XI(:,:)
  END TYPE COUPLING_PARAMETERS

  !Module variables

  TYPE (ARRAY_PROBLEM_BASE) BASE_INFO
  TYPE (ARRAY_MESH) MESH_INFO(3)
  TYPE (DARCY_PARAMETERS) DARCY
  TYPE(FIELD_TYPE), POINTER :: FIELD
  TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
  TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:)


  INTEGER(INTG), DIMENSION(:), ALLOCATABLE:: NodesPerElement
  INTEGER(INTG), DIMENSION(:,:), ALLOCATABLE::ElementNodes
  INTEGER(INTG):: NumberOfFields
  INTEGER(INTG):: NumberOfDimensions
  INTEGER(INTG):: ValueIndex
  INTEGER(INTG):: NumberOfVariableComponents
  INTEGER(INTG):: NumberOfMeshComponents
  INTEGER(INTG):: NumberOfMaterialComponents
  INTEGER(INTG):: NumberOfNodesDefined
  INTEGER(INTG):: NumberOfFieldComponent(3)
  INTEGER(INTG):: NumberOfElements
  INTEGER(INTG):: GlobalElementNumber(10)
  INTEGER(INTG):: MaxNodesPerElement
  INTEGER(INTG):: MaxNodesPerMeshComponent
  INTEGER(INTG):: ELEMENT_NUMBER
  INTEGER(INTG):: lagrange_simplex
  INTEGER(INTG), DIMENSION(:), ALLOCATABLE:: NodesPerMeshComponent
  INTEGER(INTG), DIMENSION(:), ALLOCATABLE::SimplexOutputHelp,HexOutputHelp
  INTEGER(INTG) FLD, DIMEN, OPENCMISS_INTERPOLATION(3),a,b
  INTEGER(INTG) NumberOfNodesPerElement(3), ArrayOfNodesDefined(3), NumberOfElementsDefined(3), TotalNumberOfNodes
  INTEGER(INTG), DIMENSION(:,:), ALLOCATABLE::OPENCMISS_ELEM_M,OPENCMISS_ELEM_V,OPENCMISS_ELEM_P
  INTEGER(INTG):: TRI_BASIS, TET_BASIS, QUAD_BASIS, HEX_BASIS
  INTEGER(INTG):: ALLOC_ERROR
  INTEGER(INTG):: FIELD_VAR_TYPE, var_idx

  INTEGER(INTG):: parameter_set_idx

  LOGICAL :: DN

  REAL(DP), DIMENSION(:,:), ALLOCATABLE::ElementNodesScales
  REAL(DP), DIMENSION(:), ALLOCATABLE:: XI_COORDINATES,COORDINATES
!  REAL(DP):: test
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeXValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeYValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeZValue 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeUValue 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeVValue 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeWValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeUValueORG  
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeVValueORG
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeWValueORG
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePValue2  
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeMIValue ! Mass increase for coupled elasticity Darcy INRIA model
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeMUValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeLabelValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeRHOValue  
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeA0Value
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeH0Value
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeEValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeSIGMAValue  
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeKappaValue  

  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeUValue_analytic
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeVValue_analytic 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeWValue_analytic 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePValue_analytic 

  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeUValue_error
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeVValue_error
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeWValue_error
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePValue_error

  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePerm2Value 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePerm3Value 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePerm4Value 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePerm5Value 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodePerm6Value 

  REAL(DP):: ScaleFactorsPerElementNodes(10,10)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE::OPENCMISS_NODE_COORD

  CHARACTER*2 NMs(99),KNOT
  CHARACTER*60 IN_CHAR
  CHARACTER*90 NIMZ
!   CHARACTER*30 NAMz
  CHARACTER*90 NAMz


  !Interfaces
  INTERFACE FLUID_MECHANICS_IO_READ_CMHEART
    MODULE PROCEDURE FLUID_MECHANICS_IO_READ_CMHEART1
    MODULE PROCEDURE FLUID_MECHANICS_IO_READ_CMHEART3
  END INTERFACE !FLUID_MECHANICS_IO_READ_CMHEART


  PUBLIC FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION
  PUBLIC FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS
  PUBLIC FLUID_MECHANICS_IO_READ_DATA

  PUBLIC FLUID_MECHANICS_IO_WRITE_CMGUI
  PUBLIC FLUID_MECHANICS_IO_READ_CMHEART
  PUBLIC EXPORT_CONTAINER,COUPLING_PARAMETERS
  PUBLIC FLUID_MECHANICS_IO_WRITE_ENCAS,FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS,FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS_PPE
  PUBLIC FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK


  PUBLIC FLUID_MECHANICS_IO_READ_DARCY_PARAMS
  PUBLIC DARCY

CONTAINS

  ! OK
  !================================================================================================================================
  !

  !> Writes solution into cmgui formats exelem and exnode.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_CMGUI(REGION, EQUATIONS_SET_GLOBAL_NUMBER, NAME, ERR, ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for
!     TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    CHARACTER(14) :: NAME !<the prefix name of file.
    INTEGER(INTG) :: ERR !<The error code
    INTEGER(INTG) :: EQUATIONS_SET_GLOBAL_NUMBER !<The error code
    TYPE(VARYING_STRING):: ERROR !<The error string
    !Local Variables
    INTEGER(INTG):: I,J,K,icompartment
    INTEGER(INTG):: MATERIAL_INTERPOLATION_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD !<A pointer to the equations set field
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)

    CALL ENTERS("FLUID_MECHANICS_IO_WRITE_CMGUI",ERR,ERROR,*999)

    IF (ALLOCATED(NodesPerElement)) DEALLOCATE(NodesPerElement)
    IF (ALLOCATED(NodesPerMeshComponent)) DEALLOCATE(NodesPerMeshComponent)
    IF (ALLOCATED(XI_COORDINATES)) DEALLOCATE(XI_COORDINATES)
    IF (ALLOCATED(COORDINATES)) DEALLOCATE(COORDINATES)
    IF (ALLOCATED(NodeXValue)) DEALLOCATE(NodeXValue)
    IF (ALLOCATED(NodeYValue)) DEALLOCATE(NodeYValue)
    IF (ALLOCATED(NodeZValue)) DEALLOCATE(NodeZValue)
    IF (ALLOCATED(NodeUValue)) DEALLOCATE(NodeUValue)
    IF (ALLOCATED(NodeVValue)) DEALLOCATE(NodeVValue)
    IF (ALLOCATED(NodeWValue)) DEALLOCATE(NodeWValue)
    IF (ALLOCATED(NodeUValueORG)) DEALLOCATE(NodeUValueORG)
    IF (ALLOCATED(NodeVValueORG)) DEALLOCATE(NodeVValueORG)
    IF (ALLOCATED(NodeWValueORG)) DEALLOCATE(NodeWValueORG)
    IF (ALLOCATED(NodePValue)) DEALLOCATE(NodePValue)
    IF (ALLOCATED(NodePValue2)) DEALLOCATE(NodePValue2)
    IF (ALLOCATED(NodeMIValue)) DEALLOCATE(NodeMIValue)
    IF (ALLOCATED(NodeMUValue)) DEALLOCATE(NodeMUValue)
    IF (ALLOCATED(NodeLabelValue)) DEALLOCATE(NodeLabelValue)
    IF (ALLOCATED(NodeRHOValue)) DEALLOCATE(NodeRHOValue)
    IF (ALLOCATED(NodeA0Value)) DEALLOCATE(NodeA0Value)
    IF (ALLOCATED(NodeH0Value)) DEALLOCATE(NodeH0Value)
    IF (ALLOCATED(NodeEValue)) DEALLOCATE(NodeEValue)
    IF (ALLOCATED(NodeSIGMAValue)) DEALLOCATE(NodeSIGMAValue)
    IF (ALLOCATED(NodeKappaValue)) DEALLOCATE(NodeKappaValue)
    IF (ALLOCATED(ElementNodesScales)) DEALLOCATE(ElementNodesScales)
    IF (ALLOCATED(ElementNodes)) DEALLOCATE(ElementNodes)

    IF (ALLOCATED(NodePerm2Value)) DEALLOCATE(NodePerm2Value) 
    IF (ALLOCATED(NodePerm3Value)) DEALLOCATE(NodePerm3Value) 
    IF (ALLOCATED(NodePerm4Value)) DEALLOCATE(NodePerm4Value) 
    IF (ALLOCATED(NodePerm5Value)) DEALLOCATE(NodePerm5Value) 
    IF (ALLOCATED(NodePerm6Value)) DEALLOCATE(NodePerm6Value) 

    !chrm, 20.08.09
    IF (ALLOCATED(NodeUValue_analytic)) DEALLOCATE(NodeUValue_analytic)
    IF (ALLOCATED(NodeVValue_analytic)) DEALLOCATE(NodeVValue_analytic)
    IF (ALLOCATED(NodeWValue_analytic)) DEALLOCATE(NodeWValue_analytic)
    IF (ALLOCATED(NodePValue_analytic)) DEALLOCATE(NodePValue_analytic)

    IF (ALLOCATED(NodeUValue_error)) DEALLOCATE(NodeUValue_error)
    IF (ALLOCATED(NodeVValue_error)) DEALLOCATE(NodeVValue_error)
    IF (ALLOCATED(NodeWValue_error)) DEALLOCATE(NodeWValue_error)
    IF (ALLOCATED(NodePValue_error)) DEALLOCATE(NodePValue_error)

    KNOT = '0'
    NMs(1) = '1'
    NMs(2) = '2'
    NMs(3) = '3'
    NMs(4) = '4'
    NMs(5) = '5'
    NMs(6) = '6'
    NMs(7) = '7'
    NMs(8) = '8'
    NMs(9) = '9'

    K = 9
    DO I = 1,9
      K = K + 1
      NMs(K) = TRIM(NMs(I))//TRIM(KNOT)
      DO J = 1,9
        K = K + 1
        NMs(K) = TRIM(NMs(I))//TRIM(NMs(J))
      END DO
    END DO

    EQUATIONS_SET => REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr

!---tob
!     FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE%VARIABLE_TYPE
!     ! '1' associated with linear matrix

    var_idx = 1
    FIELD_VAR_TYPE = FIELD_U_VARIABLE_TYPE
    SELECT CASE(EQUATIONS_SET%CLASS)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
          var_idx = 3
          FIELD_VAR_TYPE = FIELD_V_VARIABLE_TYPE
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
           EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
           CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
           icompartment=EQUATIONS_SET_FIELD_DATA(1)
           var_idx = 3+2*(icompartment-1)
           FIELD_VAR_TYPE=FIELD_V_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(icompartment-1))
        END SELECT
      END SELECT
    END SELECT

   


    parameter_set_idx = 1
    SELECT CASE(EQUATIONS_SET%CLASS)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
!          parameter_set_idx = 3  ! of shared dependent variable field
        END SELECT
      END SELECT
    END SELECT
!---toe

!     NumberOfFields=REGION%fields%number_of_fields
! Hack for ALE... to be removed later
    NumberOfFields=3
    NumberOfDimensions=REGION%coordinate_system%number_of_dimensions
!     NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
!       & variables(1)%number_of_components
    NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
      & variables(var_idx)%number_of_components

    NumberOfMaterialComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
      & variables(1)%number_of_components
    NumberOfElements=REGION%meshes%meshes(1)%ptr%number_of_elements
    NumberOfMeshComponents=REGION%meshes%meshes(1)%ptr%number_of_components
!     IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(NumberOfMeshComponents))

    IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(MAX(NumberOfMeshComponents,NumberOfElements)))

    IF(.NOT.ALLOCATED(NodesPerMeshComponent)) ALLOCATE(NodesPerMeshComponent(NumberOfMeshComponents))
    MaxNodesPerElement=0

    DO I=1,NumberOfMeshComponents
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(1)%basis%number_of_element_parameters
      NodesPerMeshComponent(I)=REGION%meshes%meshes(1)%ptr%topology(I)%ptr%nodes%number_of_nodes
    END DO


!     MaxNodesPerElement=NodesPerElement(1)
    MaxNodesPerMeshComponent=NodesPerMeshComponent(1)


    DO I=1,NumberOfElements
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters
      MaxNodesPerElement=MAX(NodesPerElement(1),NodesPerElement(I))
    END DO


    IF(.NOT.ALLOCATED(XI_COORDINATES))  ALLOCATE(XI_COORDINATES(NumberOfDimensions))
    IF(.NOT.ALLOCATED(COORDINATES)) ALLOCATE(COORDINATES(NumberOfDimensions))
    IF(.NOT.ALLOCATED(NodeXValue)) ALLOCATE(NodeXValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeYValue)) ALLOCATE(NodeYValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeZValue)) ALLOCATE(NodeZValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeUValue)) ALLOCATE(NodeUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue)) ALLOCATE(NodeVValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue)) ALLOCATE(NodeWValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeUValueORG)) ALLOCATE(NodeUValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValueORG)) ALLOCATE(NodeVValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValueORG)) ALLOCATE(NodeWValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue)) ALLOCATE(NodePValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue2)) ALLOCATE(NodePValue2(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeMIValue)) ALLOCATE(NodeMIValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeMUValue)) ALLOCATE(NodeMUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeLabelValue)) ALLOCATE(NodeLabelValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeRHOValue)) ALLOCATE(NodeRHOValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeA0Value)) ALLOCATE(NodeA0Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeH0Value)) ALLOCATE(NodeH0Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeEValue)) ALLOCATE(NodeEValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSIGMAValue)) ALLOCATE(NodeSIGMAValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeKappaValue)) ALLOCATE(NodeKappaValue(NodesPerMeshComponent(1)))
!     IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,NodesPerElement(1)))
    IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,MaxNodesPerElement))
!     IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,NodesPerElement(1)))
    IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,MaxNodesPerElement))

    IF(.NOT.ALLOCATED(NodePerm2Value)) ALLOCATE(NodePerm2Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm3Value)) ALLOCATE(NodePerm3Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm4Value)) ALLOCATE(NodePerm4Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm5Value)) ALLOCATE(NodePerm5Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm6Value)) ALLOCATE(NodePerm6Value(NodesPerMeshComponent(1)))

    !chrm, 20.08.09
    IF(.NOT.ALLOCATED(NodeUValue_analytic)) ALLOCATE(NodeUValue_analytic(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue_analytic)) ALLOCATE(NodeVValue_analytic(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue_analytic)) ALLOCATE(NodeWValue_analytic(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue_analytic)) ALLOCATE(NodePValue_analytic(NodesPerMeshComponent(1)))

    IF(.NOT.ALLOCATED(NodeUValue_error)) ALLOCATE(NodeUValue_error(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue_error)) ALLOCATE(NodeVValue_error(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue_error)) ALLOCATE(NodeWValue_error(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue_error)) ALLOCATE(NodePValue_error(NodesPerMeshComponent(1)))

    CALL ENTERS("CMGUI OUTPUT",ERR,ERROR,*999)

    FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field

    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATED_POINT)
    
    CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,INTERPOLATION_PARAMETERS &
    & ,ERR,ERROR,*999)
    CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)

    DO I=1,NumberOfElements

      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters

!       DO J=1,NodesPerElement(1)
      DO J=1,NodesPerElement(I)

        ELEMENT_NUMBER=I
!        XI_COORDINATES(1)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation% &
!          & geometric_interp_parameters%bases(1)%ptr%node_position_index(J,1)-1.0)/(REGION%equations_sets% &
!          & equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1) &
!          & %ptr%number_of_nodes_xi(1)-1.0)
!        XI_COORDINATES(2)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation% &
!          & geometric_interp_parameters%bases(1)%ptr%node_position_index(J,2)-1.0)/(REGION%equations_sets% &
!          & equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1) &
!          & %ptr%number_of_nodes_xi(2)-1.0)
        XI_COORDINATES(1)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,1)-1.0)/(REGION% &
          & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xic(1)-1.0)
        IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3)THEN
          XI_COORDINATES(2)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
            & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,2)-1.0)/(REGION% &
            & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
            & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xic(2)-1.0)
        END IF
        IF(NumberOfDimensions==3)THEN
          XI_COORDINATES(3)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
            & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,3)-1.0)/(REGION% &
            & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
            & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xic(3)-1.0)
        END IF

!Start: This is a hack for 3D simplex elements
        IF(NumberOfDimensions==2)THEN
          IF (NodesPerElement(1)==3) THEN
            IF(J==1)  XI_COORDINATES=(/0.0_DP,1.0_DP/)
            IF(J==2)  XI_COORDINATES=(/1.0_DP,0.0_DP/)
            IF(J==3)  XI_COORDINATES=(/1.0_DP,1.0_DP/)
          ELSE IF (NodesPerElement(1)==6) THEN
            IF(J==1)  XI_COORDINATES=(/0.0_DP,1.0_DP/)
            IF(J==2)  XI_COORDINATES=(/1.0_DP,0.0_DP/)
            IF(J==3)  XI_COORDINATES=(/1.0_DP,1.0_DP/)
            IF(J==4)  XI_COORDINATES=(/0.5_DP,0.5_DP/)
            IF(J==5)  XI_COORDINATES=(/1.0_DP,0.5_DP/)
            IF(J==6)  XI_COORDINATES=(/0.5_DP,1.0_DP/)
          ELSE IF (NodesPerElement(1)==10) THEN
            IF(J==1)  XI_COORDINATES=(/0.0_DP,1.0_DP/)
            IF(J==2)  XI_COORDINATES=(/1.0_DP,0.0_DP/)
            IF(J==3)  XI_COORDINATES=(/1.0_DP,1.0_DP/)
            IF(J==4)  XI_COORDINATES=(/1.0_DP/3.0_DP,2.0_DP/3.0_DP/)
            IF(J==5)  XI_COORDINATES=(/2.0_DP/3.0_DP,1.0_DP/3.0_DP/)
            IF(J==6)  XI_COORDINATES=(/1.0_DP,1.0_DP/3.0_DP/)
            IF(J==7)  XI_COORDINATES=(/1.0_DP,2.0_DP/3.0_DP/)
            IF(J==8)  XI_COORDINATES=(/2.0_DP/3.0_DP,1.0_DP/)
            IF(J==9)  XI_COORDINATES=(/1.0_DP/3.0_DP,1.0_DP/)
            IF(J==10)  XI_COORDINATES=(/2.0_DP/3.0_DP,2.0_DP/3.0_DP/)
          ENDIF
        ELSE IF(NumberOfDimensions==3)THEN
          IF (NodesPerElement(1)==4) THEN
            IF(J==1)  XI_COORDINATES=(/0.0_DP,1.0_DP,1.0_DP/)
            IF(J==2)  XI_COORDINATES=(/1.0_DP,0.0_DP,1.0_DP/)
            IF(J==3)  XI_COORDINATES=(/1.0_DP,1.0_DP,0.0_DP/)
            IF(J==4)  XI_COORDINATES=(/1.0_DP,1.0_DP,1.0_DP/)
          ELSE IF (NodesPerElement(1)==10) THEN
            IF(J==1)  XI_COORDINATES=(/0.0_DP,1.0_DP,1.0_DP/)
            IF(J==2)  XI_COORDINATES=(/1.0_DP,0.0_DP,1.0_DP/)
            IF(J==3)  XI_COORDINATES=(/1.0_DP,1.0_DP,0.0_DP/)
            IF(J==4)  XI_COORDINATES=(/1.0_DP,1.0_DP,1.0_DP/)
            IF(J==5)  XI_COORDINATES=(/0.5_DP,0.5_DP,1.0_DP/)
            IF(J==6)  XI_COORDINATES=(/0.5_DP,1.0_DP,0.5_DP/)
            IF(J==7)  XI_COORDINATES=(/0.5_DP,1.0_DP,1.0_DP/)
            IF(J==8)  XI_COORDINATES=(/1.0_DP,0.5_DP,0.5_DP/)
            IF(J==9)  XI_COORDINATES=(/1.0_DP,1.0_DP,0.5_DP/)
            IF(J==10)  XI_COORDINATES=(/1.0_DP,0.5_DP,1.0_DP/)
          ELSE IF (NodesPerElement(1)==20) THEN
            IF(J==1)  XI_COORDINATES=(/0.0_DP,1.0_DP,1.0_DP/)
            IF(J==2)  XI_COORDINATES=(/1.0_DP,0.0_DP,1.0_DP/)
            IF(J==3)  XI_COORDINATES=(/1.0_DP,1.0_DP,0.0_DP/)
            IF(J==4)  XI_COORDINATES=(/1.0_DP,1.0_DP,1.0_DP/)
            IF(J==5)  XI_COORDINATES=(/1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP/)
            IF(J==6)  XI_COORDINATES=(/2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP/)
            IF(J==7)  XI_COORDINATES=(/1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP/)
            IF(J==8)  XI_COORDINATES=(/2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP/)
            IF(J==9)  XI_COORDINATES=(/1.0_DP/3.0_DP,1.0_DP,1.0_DP/)
            IF(J==10)  XI_COORDINATES=(/2.0_DP/3.0_DP,1.0_DP,1.0_DP/)
            IF(J==11)  XI_COORDINATES=(/1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP/)
            IF(J==12)  XI_COORDINATES=(/1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP/)
            IF(J==13)  XI_COORDINATES=(/1.0_DP,1.0_DP,1.0_DP/3.0_DP/)
            IF(J==14)  XI_COORDINATES=(/1.0_DP,1.0_DP,2.0_DP/3.0_DP/)
            IF(J==15)  XI_COORDINATES=(/1.0_DP,1.0_DP/3.0_DP,1.0_DP/)
            IF(J==16)  XI_COORDINATES=(/1.0_DP,2.0_DP/3.0_DP,1.0_DP/)
            IF(J==17)  XI_COORDINATES=(/2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP/)
            IF(J==18)  XI_COORDINATES=(/2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP/)
            IF(J==19)  XI_COORDINATES=(/2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP/)
            IF(J==20)  XI_COORDINATES=(/1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP/)
          ENDIF
        ENDIF

!End: This is a hack for 3D simplex elements

        !K is global node number
        K=REGION%meshes%meshes(1)%ptr%topology(1)%ptr%elements%elements(I)%global_element_nodes(J)

        IF(NumberOfDimensions==3)THEN
          COORDINATES=(/1,1,1/)
        ELSE IF(NumberOfDimensions==2)THEN
          COORDINATES=(/1,1/)
        END IF

        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & INTERPOLATION_PARAMETERS(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
        NodeXValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)
!        NodeYValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1) &
!          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))
        IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3)THEN
          NodeYValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))
        END IF
        IF(NumberOfDimensions==3)THEN
          NodeZValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
        END IF

        !chrm, 19.12.2010
        !If geometry uses a quadratic mesh, but Darcy velocity and mass increase only a linear one,
        !   then we do need to interpolate the velocities and mass increase on the geometric mesh
        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE) ) THEN
          IF ((EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) .OR. &
          & (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)) THEN
            NodeUValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(1,1)
            NodeVValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(2,1)
            NodeWValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(3,1)
          ENDIF
        ELSE
!           NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field%variables(1) &
          NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
            & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)
!           NodeVValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field%variables(1) &
          IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS .OR. &
            & EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) THEN  
           NodeVValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
             & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters% &
             & cmiss%data_dp(K+NodesPerMeshComponent(1))

            IF(NumberOfDimensions==3 .OR. NumberOfDimensions==1)THEN
           NodeWValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
!                & variables(1)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
               & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters% &
               & cmiss%data_dp(K+2*NodesPerMeshComponent(1))
            END IF
          END IF
        END IF


! ! !       NodeUValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(1,1)
! ! !       NodeVValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(2,1)
! ! !       NodeWValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(3,1)


        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
          & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
          IF(NumberOfDimensions==3)THEN
            NodePValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(4,1)
          ELSE IF(NumberOfDimensions==2)THEN
            NodePValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(3,1)
          END IF
        END IF

!---tob: Mass increase for coupled elasticity Darcy INRIA model
        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE) &
            & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
              IF(NumberOfDimensions==3)THEN
                NodeMIValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(5,1)
              END IF
        END IF
!---toe

        MATERIAL_INTERPOLATION_TYPE=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
          & materials%materials_field%variables(1)%COMPONENTS(1)%INTERPOLATION_TYPE 

        IF(MATERIAL_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION)THEN
          NodeMUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)
          NodeRHOValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))

          IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS)THEN
            IF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE)THEN
              IF( (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) &
              & .OR. (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) )THEN
                NodeKappaValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
                  & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
              END IF
            END IF
          END IF

          IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
            & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) )THEN
            !--- Remaining tensor material data of permeability tensor
            NodePerm2Value(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
              & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
            NodePerm3Value(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
              & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+3*NodesPerMeshComponent(1))
            NodePerm4Value(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
              & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+4*NodesPerMeshComponent(1))
            NodePerm5Value(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
              & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+5*NodesPerMeshComponent(1))
            NodePerm6Value(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
              & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+6*NodesPerMeshComponent(1))
          END IF
        ELSE !default to FIELD_CONSTANT_INTERPOLATION
          NodeMUValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1)
          NodeRHOValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(2)

          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_1DTRANSIENT_NAVIER_STOKES_SUBTYPE) THEN
          NodeEValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(3)
          NodeH0Value=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(6)
          NodeA0Value=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(9)
          NodeSIGMAValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(12)
          END IF

          IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS)THEN
            IF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE)THEN
              IF( (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) &
              & .OR. (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) )THEN
                NodeKappaValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
                  & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(3)
              END IF
            END IF
          END IF
        END IF

      END DO 
    END DO

    DN=.TRUE.
    IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
      & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE) &
        & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) )THEN
      DN=.FALSE.
    END IF
    IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS)THEN
      DN=.FALSE.
    END IF

    IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3)THEN
      IF(DN) THEN
      ! output for DN only
        DO K=1,NodesPerMeshComponent(2)
          IF(NumberOfDimensions==2)THEN
            NodePValue2(K)=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%variables(1) &
             & %parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(2*NodesPerMeshComponent(1)+K)
          ELSE IF(NumberOfDimensions==3)THEN
            NodePValue2(K)=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%variables(1) &
              & %parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(3*NodesPerMeshComponent(1)+K)
          ENDIF
        ENDDO
      ENDIF
    ENDIF

    IF( NumberOfDimensions==3 )THEN
      !For 3D, the following call works ...
      lagrange_simplex=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations% &
        & interpolation%geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%type
!         lagrange_simplex=2

    ELSE IF (NumberOfDimensions==1) THEN
      lagrange_simplex=1
    ELSE IF (NumberOfDimensions==2) THEN
      !chrm, 20.08.09:
      ! ... but the above call does not work for 2D.
      !Thus, for 2D, we hard-wire it to 'quad':
      IF(MaxNodesPerElement==4.OR.MaxNodesPerElement==9.OR.MaxNodesPerElement==16) THEN
        lagrange_simplex=1
      ELSE IF(MaxNodesPerElement==3.OR.MaxNodesPerElement==6.OR.MaxNodesPerElement==10) THEN
        lagrange_simplex=2
      ENDIF
    END IF

    NumberOfFieldComponent(1)=NumberOfDimensions
    NumberOfFieldComponent(2)=NumberOfVariableComponents
    NumberOfFieldComponent(3)=NumberOfMaterialComponents

    DO I=1,NumberOfElements
!       DO J=1,NodesPerElement(1)
      DO J=1,NodesPerElement(I)
        ElementNodes(I,J)=REGION%meshes%meshes(1)%ptr%topology(1)% &
          & ptr%elements%elements(I)%global_element_nodes(J)
        ElementNodesScales(I,J)=1.0000000000000000E+00
      END DO
    END DO

    CALL FLUID_MECHANICS_IO_WRITE_NODES_CMGUI(NAME,EQUATIONS_SET)
    CALL FLUID_MECHANICS_IO_WRITE_ELEMENTS_CMGUI(NAME)

    CALL EXITS("FLUID_MECHANICS_IO_WRITE_CMGUI")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_CMGUI",ERR,ERROR)    
    CALL EXITS("FLUID_MECHANICS_IO_WRITE_CMGUI")
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_CMGUI


  ! OK
  !================================================================================================================================
  !  

  !> Writes solution into encas
!   SUBROUTINE FLUID_MECHANICS_IO_WRITE_ENCAS(REGION, NAME, ERR, ERROR,*)
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_ENCAS(REGION, EQUATIONS_SET_GLOBAL_NUMBER, NAME, ERR, ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for
!     TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    CHARACTER(14) :: NAME !<the prefix name of file.
    INTEGER(INTG) :: EQUATIONS_SET_GLOBAL_NUMBER !<The error code
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING):: ERROR !<The error string
    !Local Variables
    INTEGER(INTG):: I,J,K,icompartment
    INTEGER(INTG):: MATERIAL_INTERPOLATION_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD !<A pointer to the equations set field
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)

    CALL ENTERS("FLUID_MECHANICS_IO_WRITE_ENCAS",ERR,ERROR,*999)

    IF (ALLOCATED(NodesPerElement)) DEALLOCATE(NodesPerElement)
    IF (ALLOCATED(NodesPerMeshComponent)) DEALLOCATE(NodesPerMeshComponent)
    IF (ALLOCATED(XI_COORDINATES)) DEALLOCATE(XI_COORDINATES)
    IF (ALLOCATED(COORDINATES)) DEALLOCATE(COORDINATES)
    IF (ALLOCATED(NodeXValue)) DEALLOCATE(NodeXValue)
    IF (ALLOCATED(NodeYValue)) DEALLOCATE(NodeYValue)
    IF (ALLOCATED(NodeZValue)) DEALLOCATE(NodeZValue)
    IF (ALLOCATED(NodeUValue)) DEALLOCATE(NodeUValue)
    IF (ALLOCATED(NodeVValue)) DEALLOCATE(NodeVValue)
    IF (ALLOCATED(NodeWValue)) DEALLOCATE(NodeWValue)
    IF (ALLOCATED(NodeUValueORG)) DEALLOCATE(NodeUValueORG)
    IF (ALLOCATED(NodeVValueORG)) DEALLOCATE(NodeVValueORG)
    IF (ALLOCATED(NodeWValueORG)) DEALLOCATE(NodeWValueORG)
    IF (ALLOCATED(NodePValue)) DEALLOCATE(NodePValue)
    IF (ALLOCATED(NodePValue2)) DEALLOCATE(NodePValue2)
    IF (ALLOCATED(NodeMIValue)) DEALLOCATE(NodeMIValue)
    IF (ALLOCATED(NodeMUValue)) DEALLOCATE(NodeMUValue)
    IF (ALLOCATED(NodeLabelValue)) DEALLOCATE(NodeLabelValue)
    IF (ALLOCATED(NodeRHOValue)) DEALLOCATE(NodeRHOValue)
    IF (ALLOCATED(NodeA0Value)) DEALLOCATE(NodeA0Value)
    IF (ALLOCATED(NodeH0Value)) DEALLOCATE(NodeH0Value)
    IF (ALLOCATED(NodeEValue)) DEALLOCATE(NodeEValue)
    IF (ALLOCATED(NodeSIGMAValue)) DEALLOCATE(NodeSIGMAValue)
    IF (ALLOCATED(NodeKappaValue)) DEALLOCATE(NodeKappaValue)
    IF (ALLOCATED(ElementNodesScales)) DEALLOCATE(ElementNodesScales)
    IF (ALLOCATED(ElementNodes)) DEALLOCATE(ElementNodes)

    !chrm, 20.08.09
    IF (ALLOCATED(NodeUValue_analytic)) DEALLOCATE(NodeUValue_analytic)
    IF (ALLOCATED(NodeVValue_analytic)) DEALLOCATE(NodeVValue_analytic)
    IF (ALLOCATED(NodeWValue_analytic)) DEALLOCATE(NodeWValue_analytic)
    IF (ALLOCATED(NodePValue_analytic)) DEALLOCATE(NodePValue_analytic)

    IF (ALLOCATED(NodeUValue_error)) DEALLOCATE(NodeUValue_error)
    IF (ALLOCATED(NodeVValue_error)) DEALLOCATE(NodeVValue_error)
    IF (ALLOCATED(NodeWValue_error)) DEALLOCATE(NodeWValue_error)
    IF (ALLOCATED(NodePValue_error)) DEALLOCATE(NodePValue_error)

    KNOT = '0'
    NMs(1) = '1'
    NMs(2) = '2'
    NMs(3) = '3'
    NMs(4) = '4'
    NMs(5) = '5'
    NMs(6) = '6'
    NMs(7) = '7'
    NMs(8) = '8'
    NMs(9) = '9'

    K = 9
    DO I = 1,9
      K = K + 1
      NMs(K) = TRIM(NMs(I))//TRIM(KNOT)
      DO J = 1,9
        K = K + 1
        NMs(K) = TRIM(NMs(I))//TRIM(NMs(J))
      END DO
    END DO

    EQUATIONS_SET => REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr

!---tob
!     FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE%VARIABLE_TYPE
!     ! '1' associated with linear matrix

    var_idx = 1
    FIELD_VAR_TYPE = FIELD_U_VARIABLE_TYPE
    SELECT CASE(EQUATIONS_SET%CLASS)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
          var_idx = 3
          FIELD_VAR_TYPE = FIELD_V_VARIABLE_TYPE
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
           EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
           CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
           icompartment=EQUATIONS_SET_FIELD_DATA(1)
           var_idx = 3+2*(icompartment-1)
           FIELD_VAR_TYPE=FIELD_V_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(icompartment-1))
        END SELECT
      END SELECT
    END SELECT

    parameter_set_idx = 1 
    SELECT CASE(EQUATIONS_SET%CLASS)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
!          parameter_set_idx = 3  ! of shared dependent variable field
        END SELECT
      END SELECT
    END SELECT
!---toe

!     NumberOfFields=REGION%fields%number_of_fields
! Hack for ALE... to be removed later
    NumberOfFields=3
    NumberOfDimensions=REGION%coordinate_system%number_of_dimensions
!     NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
!       & variables(1)%number_of_components
    NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
      & variables(var_idx)%number_of_components

    NumberOfMaterialComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
      & variables(1)%number_of_components
    NumberOfElements=REGION%meshes%meshes(1)%ptr%number_of_elements
    NumberOfMeshComponents=REGION%meshes%meshes(1)%ptr%number_of_components
!     IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(NumberOfMeshComponents))

    IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(MAX(NumberOfMeshComponents,NumberOfElements)))

    IF(.NOT.ALLOCATED(NodesPerMeshComponent)) ALLOCATE(NodesPerMeshComponent(NumberOfMeshComponents))
    MaxNodesPerElement=0

    DO I=1,NumberOfMeshComponents
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(1)%basis%number_of_element_parameters
      NodesPerMeshComponent(I)=REGION%meshes%meshes(1)%ptr%topology(I)%ptr%nodes%number_of_nodes
    END DO


!     MaxNodesPerElement=NodesPerElement(1)
    MaxNodesPerMeshComponent=NodesPerMeshComponent(1)


    DO I=1,NumberOfElements
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters
      MaxNodesPerElement=MAX(NodesPerElement(1),NodesPerElement(I))
    END DO


    IF(.NOT.ALLOCATED(XI_COORDINATES))  ALLOCATE(XI_COORDINATES(NumberOfDimensions))
    IF(.NOT.ALLOCATED(COORDINATES)) ALLOCATE(COORDINATES(NumberOfDimensions))
    IF(.NOT.ALLOCATED(NodeXValue)) ALLOCATE(NodeXValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeYValue)) ALLOCATE(NodeYValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeZValue)) ALLOCATE(NodeZValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeUValue)) ALLOCATE(NodeUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue)) ALLOCATE(NodeVValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue)) ALLOCATE(NodeWValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeUValueORG)) ALLOCATE(NodeUValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValueORG)) ALLOCATE(NodeVValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValueORG)) ALLOCATE(NodeWValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue)) ALLOCATE(NodePValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue2)) ALLOCATE(NodePValue2(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeMIValue)) ALLOCATE(NodeMIValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeMUValue)) ALLOCATE(NodeMUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeLabelValue)) ALLOCATE(NodeLabelValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeRHOValue)) ALLOCATE(NodeRHOValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeA0Value)) ALLOCATE(NodeA0Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeH0Value)) ALLOCATE(NodeH0Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeEValue)) ALLOCATE(NodeEValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSIGMAValue)) ALLOCATE(NodeSIGMAValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeKappaValue)) ALLOCATE(NodeKappaValue(NodesPerMeshComponent(1)))
!     IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,NodesPerElement(1)))
    IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,MaxNodesPerElement))
!     IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,NodesPerElement(1)))
    IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,MaxNodesPerElement))

    !chrm, 20.08.09
    IF(.NOT.ALLOCATED(NodeUValue_analytic)) ALLOCATE(NodeUValue_analytic(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue_analytic)) ALLOCATE(NodeVValue_analytic(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue_analytic)) ALLOCATE(NodeWValue_analytic(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue_analytic)) ALLOCATE(NodePValue_analytic(NodesPerMeshComponent(1)))

    IF(.NOT.ALLOCATED(NodeUValue_error)) ALLOCATE(NodeUValue_error(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue_error)) ALLOCATE(NodeVValue_error(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue_error)) ALLOCATE(NodeWValue_error(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue_error)) ALLOCATE(NodePValue_error(NodesPerMeshComponent(1)))

    CALL ENTERS("CMGUI OUTPUT",ERR,ERROR,*999)

    FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field

    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATED_POINT)
    
    CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,INTERPOLATION_PARAMETERS &
    & ,ERR,ERROR,*999)
    CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)

    DO I=1,NumberOfElements

      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters

!       DO J=1,NodesPerElement(1)
      DO J=1,NodesPerElement(I)

        ELEMENT_NUMBER=I
        XI_COORDINATES(1)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,1)-1.0)/(REGION% &
          & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xic(1)-1.0)
        XI_COORDINATES(2)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,2)-1.0)/(REGION% &
          & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xic(2)-1.0)
        XI_COORDINATES(3)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,3)-1.0)/(REGION% &
          & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
          & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xic(3)-1.0)
        IF(NumberOfDimensions==2)THEN
          STOP 'Encas format only available for 3D hex and tets'
        END IF

        !K is global node number
        K=REGION%meshes%meshes(1)%ptr%topology(1)%ptr%elements%elements(I)%global_element_nodes(J)

        COORDINATES=(/1,1,1/)

        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & INTERPOLATION_PARAMETERS(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
        NodeXValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)
        NodeYValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))

        IF(NumberOfDimensions==3)THEN
          NodeZValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
        END IF

!         NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field%variables(1) &
        NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
          & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)
!         NodeVValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field%variables(1) &
        NodeVValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
          & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters% &
          & cmiss%data_dp(K+NodesPerMeshComponent(1))

        IF(NumberOfDimensions==3)THEN
          NodeWValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
!             & variables(1)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
            & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters% &
            & cmiss%data_dp(K+2*NodesPerMeshComponent(1))
        END IF

! ! !       NodeUValue(K)=INTERPOLATED_POINT%VALUES(1,1)
! ! !       NodeVValue(K)=INTERPOLATED_POINT%VALUES(2,1)
! ! !       NodeWValue(K)=INTERPOLATED_POINT%VALUES(3,1)

        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) &
          & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
          IF(NumberOfDimensions==3)THEN
            NodePValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(4,1)
          ELSE IF(NumberOfDimensions==2)THEN
            NodePValue(K)=INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr%VALUES(3,1)
          END IF
        END IF

        MATERIAL_INTERPOLATION_TYPE=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
          & materials%materials_field%variables(1)%COMPONENTS(1)%INTERPOLATION_TYPE 

        IF(MATERIAL_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION)THEN
          NodeMUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)
          NodeRHOValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))

          IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS)THEN
            IF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE)THEN
              IF( (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) &
              & .OR. (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) )THEN
                NodeKappaValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
                  & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
              END IF
            END IF
          END IF
        ELSE !default to FIELD_CONSTANT_INTERPOLATION
          NodeMUValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1)
          NodeRHOValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(2)
          IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS)THEN
            IF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE)THEN
              IF( (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) )THEN
                NodeKappaValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
                  & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(3)
              END IF
            END IF
          END IF
        END IF

      END DO 
    END DO

    IF( NumberOfDimensions==3 )THEN
      !For 3D, the following call works ...
      lagrange_simplex=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations% &
        & interpolation%geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%type
!         lagrange_simplex=2
    ELSE IF ( NumberOfDimensions==1 )THEN
      lagrange_simplex=1
    ELSE IF ( NumberOfDimensions==2 )THEN
      !chrm, 20.08.09:
      ! ... but the above call does not work for 2D.
      !Thus, for 2D, we hard-wire it to 'quad':
      IF(MaxNodesPerElement==4.OR.MaxNodesPerElement==9.OR.MaxNodesPerElement==16) THEN
        lagrange_simplex=1
      ELSE IF(MaxNodesPerElement==3.OR.MaxNodesPerElement==6.OR.MaxNodesPerElement==10) THEN
        lagrange_simplex=2
      ENDIF
    END IF

    NumberOfFieldComponent(1)=NumberOfDimensions
    NumberOfFieldComponent(2)=NumberOfVariableComponents
    NumberOfFieldComponent(3)=NumberOfMaterialComponents

    DO I=1,NumberOfElements
!       DO J=1,NodesPerElement(1)
      DO J=1,NodesPerElement(I)
        ElementNodes(I,J)=REGION%meshes%meshes(1)%ptr%topology(1)% &
          & ptr%elements%elements(I)%global_element_nodes(J)
        ElementNodesScales(I,J)=1.0000000000000000E+00
      END DO
    END DO

    CALL FLUID_MECHANICS_IO_WRITE_DATA_ENCAS(NAME,EQUATIONS_SET)


    CALL EXITS("FLUID_MECHANICS_IO_WRITE_ENCAS")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_ENCAS",ERR,ERROR)    
    CALL EXITS("FLUID_MECHANICS_IO_WRITE_ENCAS")
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_ENCAS



  ! OK
  !================================================================================================================================
  !  

  !> Executes nodes writing process.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_DATA_ENCAS(NAME,EQUATIONS_SET)

    IMPLICIT NONE
    CHARACTER(14), INTENT(IN) :: NAME !<the prefix name of file.
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
!     CHARACTER :: FILENAME !<the prefix name of file.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER:: I,K
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR 
    DOUBLE PRECISION:: velocity_magnitude
!     CHARACTER(80):: OUTPUT_FILE

  !==================
    FILENAME="./output/"//NAME//".geo"
    OPEN(UNIT=80, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(80,*)'OpenCMISS Exported Encas Model Geometry File' 
    WRITE(80,*)'' 
    WRITE(80,*)'node id assign'
    WRITE(80,*)'element id assign'
    WRITE(80,*)'part'
    WRITE(80,*)'        1'
    WRITE(80,*)'OpenCMISS'
    WRITE(80,*)'coordinates'
    WRITE(80,'(i10)')NodesPerMeshComponent(1)

    DO I = 1,NodesPerMeshComponent(1)
      WRITE(80,'(e12.5)')NodeXValue(I)
    ENDDO
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(80,'(e12.5)')NodeYValue(I)
    ENDDO
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(80,'(e12.5)')NodeZValue(I)
    ENDDO
    IF(lagrange_simplex==2) THEN !simplex
      IF(.NOT.ALLOCATED(SimplexOutputHelp)) ALLOCATE(SimplexOutputHelp(NodesPerElement(1)))
      IF(MaxNodesPerElement==4) THEN
        WRITE(80,*)'tetra4'
        WRITE(80,'(i10)')NumberOfElements
        DO K = 1,NumberOfElements
          SimplexOutputHelp(1)=ElementNodes(K,1)
          SimplexOutputHelp(2)=ElementNodes(K,2)
          SimplexOutputHelp(3)=ElementNodes(K,4)
          SimplexOutputHelp(4)=ElementNodes(K,3)
          WRITE(80,'(4(i10))')SimplexOutputHelp
        END DO
      ELSE IF(MaxNodesPerElement==10) THEN
        WRITE(80,*)'tetra10'
        WRITE(80,'(i10)')NumberOfElements
        DO K = 1,NumberOfElements
          SimplexOutputHelp(1)=ElementNodes(K,1)
          SimplexOutputHelp(2)=ElementNodes(K,2)
          SimplexOutputHelp(3)=ElementNodes(K,4)
          SimplexOutputHelp(4)=ElementNodes(K,3)
          SimplexOutputHelp(5)=ElementNodes(K,5)
          SimplexOutputHelp(6)=ElementNodes(K,10)
          SimplexOutputHelp(7)=ElementNodes(K,7)
          SimplexOutputHelp(8)=ElementNodes(K,6)
          SimplexOutputHelp(9)=ElementNodes(K,8)
          SimplexOutputHelp(10)=ElementNodes(K,9)
          WRITE(80,'(10(i10))')SimplexOutputHelp
        END DO
      ELSE
        STOP 'Encas format only available for 3D hex and tets'
      ENDIF
    ELSE IF (lagrange_simplex==1) THEN !hexa
      IF(.NOT.ALLOCATED(HexOutputHelp)) ALLOCATE(HexOutputHelp(NodesPerElement(1)))
      IF(MaxNodesPerElement==8) THEN
        WRITE(80,*)'hexa8'
        WRITE(80,'(i10)')NumberOfElements
        DO K = 1,NumberOfElements
          HexOutputHelp(1)=ElementNodes(K,1)
          HexOutputHelp(2)=ElementNodes(K,2)
          HexOutputHelp(3)=ElementNodes(K,4)
          HexOutputHelp(4)=ElementNodes(K,3)
          HexOutputHelp(5)=ElementNodes(K,5)
          HexOutputHelp(6)=ElementNodes(K,6)
          HexOutputHelp(7)=ElementNodes(K,8)
          HexOutputHelp(8)=ElementNodes(K,7)
          WRITE(80,'(8(i10))')HexOutputHelp
        END DO
      ELSE IF(MaxNodesPerElement==27) THEN
        WRITE(80,*)'hexa8'
        WRITE(80,'(i10)')8*NumberOfElements
        DO K = 1,NumberOfElements
          HexOutputHelp(1)=ElementNodes(K,1)
          HexOutputHelp(2)=ElementNodes(K,2)
          HexOutputHelp(3)=ElementNodes(K,5)
          HexOutputHelp(4)=ElementNodes(K,4)
          HexOutputHelp(5)=ElementNodes(K,10)
          HexOutputHelp(6)=ElementNodes(K,11)
          HexOutputHelp(7)=ElementNodes(K,14)
          HexOutputHelp(8)=ElementNodes(K,13)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,2)
          HexOutputHelp(2)=ElementNodes(K,3)
          HexOutputHelp(3)=ElementNodes(K,6)
          HexOutputHelp(4)=ElementNodes(K,5)
          HexOutputHelp(5)=ElementNodes(K,11)
          HexOutputHelp(6)=ElementNodes(K,12)
          HexOutputHelp(7)=ElementNodes(K,15)
          HexOutputHelp(8)=ElementNodes(K,14)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,4)
          HexOutputHelp(2)=ElementNodes(K,5)
          HexOutputHelp(3)=ElementNodes(K,8)
          HexOutputHelp(4)=ElementNodes(K,7)
          HexOutputHelp(5)=ElementNodes(K,13)
          HexOutputHelp(6)=ElementNodes(K,14)
          HexOutputHelp(7)=ElementNodes(K,17)
          HexOutputHelp(8)=ElementNodes(K,16)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,5)
          HexOutputHelp(2)=ElementNodes(K,6)
          HexOutputHelp(3)=ElementNodes(K,9)
          HexOutputHelp(4)=ElementNodes(K,8)
          HexOutputHelp(5)=ElementNodes(K,14)
          HexOutputHelp(6)=ElementNodes(K,15)
          HexOutputHelp(7)=ElementNodes(K,18)
          HexOutputHelp(8)=ElementNodes(K,17)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,10)
          HexOutputHelp(2)=ElementNodes(K,11)
          HexOutputHelp(3)=ElementNodes(K,14)
          HexOutputHelp(4)=ElementNodes(K,13)
          HexOutputHelp(5)=ElementNodes(K,19)
          HexOutputHelp(6)=ElementNodes(K,20)
          HexOutputHelp(7)=ElementNodes(K,23)
          HexOutputHelp(8)=ElementNodes(K,22)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,11)
          HexOutputHelp(2)=ElementNodes(K,12)
          HexOutputHelp(3)=ElementNodes(K,15)
          HexOutputHelp(4)=ElementNodes(K,14)
          HexOutputHelp(5)=ElementNodes(K,20)
          HexOutputHelp(6)=ElementNodes(K,21)
          HexOutputHelp(7)=ElementNodes(K,24)
          HexOutputHelp(8)=ElementNodes(K,23)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,13)
          HexOutputHelp(2)=ElementNodes(K,14)
          HexOutputHelp(3)=ElementNodes(K,17)
          HexOutputHelp(4)=ElementNodes(K,16)
          HexOutputHelp(5)=ElementNodes(K,22)
          HexOutputHelp(6)=ElementNodes(K,23)
          HexOutputHelp(7)=ElementNodes(K,26)
          HexOutputHelp(8)=ElementNodes(K,25)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
          HexOutputHelp(1)=ElementNodes(K,14)
          HexOutputHelp(2)=ElementNodes(K,15)
          HexOutputHelp(3)=ElementNodes(K,18)
          HexOutputHelp(4)=ElementNodes(K,17)
          HexOutputHelp(5)=ElementNodes(K,23)
          HexOutputHelp(6)=ElementNodes(K,24)
          HexOutputHelp(7)=ElementNodes(K,27)
          HexOutputHelp(8)=ElementNodes(K,26)
          WRITE(80,'(8(i10))')HexOutputHelp(1:8)
        END DO
      ELSE
        STOP 'Encas format only available for 3D hex and tets'
      ENDIF
    ENDIF
    CLOSE(80)
  !==================
    IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) &
      & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
            & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
            & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
      FILENAME="./output/"//NAME//".scl1"
      OPEN(UNIT=81, FILE=CHAR(FILENAME),STATUS='unknown')
      WRITE(81,*)'Absolute Pressure' 
      WRITE(81,*)'part'
      WRITE(81,*)'        1'
      WRITE(81,*)'coordinates'
      DO I = 1,NodesPerMeshComponent(1)
        WRITE(81,'(e12.5)')NodePValue(I)
      ENDDO
      CLOSE(81)
    ENDIF
 
  !==================
    FILENAME="./output/"//NAME//".scl2"
    OPEN(UNIT=82, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(82,*)'Velocity Magnitude' 
    WRITE(82,*)'part'
    WRITE(82,*)'        1'
    WRITE(82,*)'coordinates'
    DO I = 1,NodesPerMeshComponent(1)
      velocity_magnitude=sqrt(NodeUValue(I)*NodeUValue(I)+ & 
        & NodeVValue(I)*NodeVValue(I)+NodeWValue(I)*NodeWValue(I))
      WRITE(82,'(e12.5)')velocity_magnitude
    ENDDO
    CLOSE(82)
  
  !==================
    FILENAME="./output/"//NAME//".scl3"
    OPEN(UNIT=83, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(83,*)'X Velocity' 
    WRITE(83,*)'part'
    WRITE(83,*)'        1'
    WRITE(83,*)'coordinates'
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(83,'(e12.5)')NodeUValue(I)
    ENDDO
    CLOSE(83)
  
  !==================
    FILENAME="./output/"//NAME//".scl4"
    OPEN(UNIT=84, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(84,*)'Y Velocity' 
    WRITE(84,*)'part'
    WRITE(84,*)'        1'
    WRITE(84,*)'coordinates'
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(84,'(e12.5)')NodeVValue(I)
    ENDDO
    CLOSE(84)

  !==================
    FILENAME="./output/"//NAME//".scl5"
    OPEN(UNIT=85, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(85,*)'Z Velocity' 
    WRITE(85,*)'part'
    WRITE(85,*)'        1'
    WRITE(85,*)'coordinates'
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(85,'(e12.5)')NodeWValue(I)
    ENDDO
    CLOSE(85)
  !==================
    FILENAME="./output/"//NAME//".vel"
    OPEN(UNIT=86, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(86,*)'Velocity' 
    WRITE(86,*)'part'
    WRITE(86,*)'        1'
    WRITE(86,*)'coordinates'
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(86,'(e12.5)')NodeUValue(I)
    ENDDO
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(86,'(e12.5)')NodeVValue(I)
    ENDDO
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(86,'(e12.5)')NodeWValue(I)
    ENDDO
    CLOSE(86)
  !==================
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Encas data...",ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_DATA_ENCAS",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_DATA_ENCAS

  ! OK
  !================================================================================================================================
  !


  !> Executes nodes writing process.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS_PPE(NAME,start_time_step,number_of_timesteps,time_increment)

    IMPLICIT NONE

    INTEGER:: I,J,start_time_step,number_of_timesteps
    DOUBLE PRECISION:: time_increment,time

    CHARACTER(9), INTENT(IN) :: NAME !<the prefix name of file.
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
!     CHARACTER :: FILENAME !<the prefix name of file.
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR 


    FILENAME="./output/"//NAME//".case"
    OPEN(UNIT=87, FILE=CHAR(FILENAME),STATUS='unknown')

    WRITE(87,*)'FORMAT' 
    WRITE(87,*)'type:  ensight gold'
    WRITE(87,*)''
    WRITE(87,*)'GEOMETRY'
    WRITE(87,*)'model:      ./OpenCMISS.geo'
    WRITE(87,*)''
    WRITE(87,*)'VARIABLE'
    WRITE(87,*)''
    WRITE(87,*)'scalar per node:  1  Magnitude OpenCMISS.Magnitude_**'
    WRITE(87,*)'scalar per node:  1  Speed_SumSquares OpenCMISS.Speed_SumSquares_**'
    WRITE(87,*)'scalar per node:  1  Pressure OpenCMISS.Pressure_**'
    WRITE(87,*)'vector per node:  1  Velocity OpenCMISS.Velocity_**'  
    WRITE(87,*)''
    WRITE(87,*)'TIME'
    WRITE(87,*)'time set: 1 Model'
    WRITE(87,*)'number of steps: ',  number_of_timesteps
    WRITE(87,*)'filename start number:      ',  start_time_step
    WRITE(87,*)'filename increment:         1'
    WRITE(87,'(" time values:")',ADVANCE="NO")
    J=1
    DO I=1,number_of_timesteps
      J=J+1
      time=I*time_increment
      WRITE(87,'(e13.5)',ADVANCE="NO")time
      IF(J==8) THEN
        WRITE(87,*) ' '
        J=0
      ENDIF
    ENDDO
  
  CLOSE(87)
  

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Encas PPE master...",ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS_PPE",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS_PPE


  ! OK
  !================================================================================================================================
  !


  !> Executes nodes writing process.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS(NAME,start_time_step,number_of_timesteps,time_increment)

    IMPLICIT NONE

    INTEGER:: I,J,start_time_step,number_of_timesteps
    DOUBLE PRECISION:: time_increment,time

    CHARACTER(14), INTENT(IN) :: NAME !<the prefix name of file.
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
!     CHARACTER :: FILENAME !<the prefix name of file.
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR 


    FILENAME="./output/"//NAME//".encas"
    OPEN(UNIT=87, FILE=CHAR(FILENAME),STATUS='unknown')

    WRITE(87,*)'FORMAT' 
    WRITE(87,*)'type:  ensight gold'
    WRITE(87,*)''
    WRITE(87,*)'GEOMETRY'
    WRITE(87,*)'model:  1   ./TIME_STEP_****.geo'
    WRITE(87,*)''
    WRITE(87,*)'VARIABLE'
    WRITE(87,*)'scalar per node:  1  pressure                  ./TIME_STEP_****.scl1'
    WRITE(87,*)'scalar per node:  1  velocity-magnitude        ./TIME_STEP_****.scl2'
    WRITE(87,*)'scalar per node:  1  velocity-u                ./TIME_STEP_****.scl3'
    WRITE(87,*)'scalar per node:  1  velocity-v                ./TIME_STEP_****.scl4'
    WRITE(87,*)'scalar per node:  1  velocity-w                ./TIME_STEP_****.scl5'
    WRITE(87,*)'vector per node:  1  velocity                  ./TIME_STEP_****.vel'  
    WRITE(87,*)''
    WRITE(87,*)'TIME'
    WRITE(87,*)'time set: 1 Model'
    WRITE(87,*)'number of steps: ',  number_of_timesteps
    WRITE(87,*)'filename start number:      ',  start_time_step
    WRITE(87,*)'filename increment:         1'
    WRITE(87,'(" time values:")',ADVANCE="NO")
    J=1
    DO I=1,number_of_timesteps
      J=J+1
      time=I*time_increment
      WRITE(87,'(e13.5)',ADVANCE="NO")time
      IF(J==8) THEN
        WRITE(87,*) ' '
        J=0
      ENDIF
    ENDDO
  
  CLOSE(87)
  

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Encas master...",ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_MASTER_ENCAS

  ! OK
  !================================================================================================================================
  !  

  !> Writes solution into encas --- BUT FOR PRESSURE POISSON ONLY --- not for general use!
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK(REGION, EQUATIONS_SET_GLOBAL_NUMBER, NAME, ERR, ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for
!     TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    CHARACTER(14) :: NAME !<the prefix name of file.
    INTEGER(INTG) :: EQUATIONS_SET_GLOBAL_NUMBER !<The error code
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING):: ERROR !<The error string
    !Local Variables
    INTEGER(INTG):: I,J,K,icompartment
    INTEGER(INTG):: MATERIAL_INTERPOLATION_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD !<A pointer to the equations set field
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)

    CALL ENTERS("FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK",ERR,ERROR,*999)

    IF (ALLOCATED(NodesPerElement)) DEALLOCATE(NodesPerElement)
    IF (ALLOCATED(NodesPerMeshComponent)) DEALLOCATE(NodesPerMeshComponent)
    IF (ALLOCATED(XI_COORDINATES)) DEALLOCATE(XI_COORDINATES)
    IF (ALLOCATED(COORDINATES)) DEALLOCATE(COORDINATES)
    IF (ALLOCATED(NodeXValue)) DEALLOCATE(NodeXValue)
    IF (ALLOCATED(NodeYValue)) DEALLOCATE(NodeYValue)
    IF (ALLOCATED(NodeZValue)) DEALLOCATE(NodeZValue)
    IF (ALLOCATED(NodeUValue)) DEALLOCATE(NodeUValue)
    IF (ALLOCATED(NodeVValue)) DEALLOCATE(NodeVValue)
    IF (ALLOCATED(NodeWValue)) DEALLOCATE(NodeWValue)
    IF (ALLOCATED(NodeUValueORG)) DEALLOCATE(NodeUValueORG)
    IF (ALLOCATED(NodeWValueORG)) DEALLOCATE(NodeVValueORG)
    IF (ALLOCATED(NodeVValueORG)) DEALLOCATE(NodeWValueORG)
    IF (ALLOCATED(NodePValue)) DEALLOCATE(NodePValue)
    IF (ALLOCATED(NodeMUValue)) DEALLOCATE(NodeMUValue)
    IF (ALLOCATED(NodeLabelValue)) DEALLOCATE(NodeLabelValue)
    IF (ALLOCATED(NodeRHOValue)) DEALLOCATE(NodeRHOValue)
    IF (ALLOCATED(NodeKappaValue)) DEALLOCATE(NodeKappaValue)
    IF (ALLOCATED(ElementNodesScales)) DEALLOCATE(ElementNodesScales)
    IF (ALLOCATED(ElementNodes)) DEALLOCATE(ElementNodes)

    KNOT = '0'
    NMs(1) = '1'
    NMs(2) = '2'
    NMs(3) = '3'
    NMs(4) = '4'
    NMs(5) = '5'
    NMs(6) = '6'
    NMs(7) = '7'
    NMs(8) = '8'
    NMs(9) = '9'

    K = 9
    DO I = 1,9
      K = K + 1
      NMs(K) = TRIM(NMs(I))//TRIM(KNOT)
      DO J = 1,9
        K = K + 1
        NMs(K) = TRIM(NMs(I))//TRIM(NMs(J))
      END DO
    END DO

    EQUATIONS_SET => REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr

!---tob
!     FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE%VARIABLE_TYPE
!     ! '1' associated with linear matrix

    var_idx = 1
    FIELD_VAR_TYPE = FIELD_U_VARIABLE_TYPE
    SELECT CASE(EQUATIONS_SET%CLASS)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
          var_idx = 3
          FIELD_VAR_TYPE = FIELD_V_VARIABLE_TYPE
        END SELECT
      END SELECT
    END SELECT

    parameter_set_idx = 1 
    SELECT CASE(EQUATIONS_SET%CLASS)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
!          parameter_set_idx = 3  ! of shared dependent variable field
        END SELECT
      END SELECT
    END SELECT
!---toe

!     NumberOfFields=REGION%fields%number_of_fields
! Hack for ALE... to be removed later
    NumberOfFields=3
    NumberOfDimensions=REGION%coordinate_system%number_of_dimensions
!     NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
!       & variables(1)%number_of_components
    NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
      & variables(var_idx)%number_of_components

    NumberOfMaterialComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
      & variables(1)%number_of_components
    NumberOfElements=REGION%meshes%meshes(1)%ptr%number_of_elements
    NumberOfMeshComponents=REGION%meshes%meshes(1)%ptr%number_of_components
!     IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(NumberOfMeshComponents))

    IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(MAX(NumberOfMeshComponents,NumberOfElements)))

    IF(.NOT.ALLOCATED(NodesPerMeshComponent)) ALLOCATE(NodesPerMeshComponent(NumberOfMeshComponents))
    MaxNodesPerElement=0

    DO I=1,NumberOfMeshComponents
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(1)%basis%number_of_element_parameters
      NodesPerMeshComponent(I)=REGION%meshes%meshes(1)%ptr%topology(I)%ptr%nodes%number_of_nodes
    END DO


!     MaxNodesPerElement=NodesPerElement(1)
    MaxNodesPerMeshComponent=NodesPerMeshComponent(1)


    DO I=1,NumberOfElements
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters
      MaxNodesPerElement=MAX(NodesPerElement(1),NodesPerElement(I))
    END DO


    IF(.NOT.ALLOCATED(XI_COORDINATES))  ALLOCATE(XI_COORDINATES(NumberOfDimensions))
    IF(.NOT.ALLOCATED(COORDINATES)) ALLOCATE(COORDINATES(NumberOfDimensions))
    IF(.NOT.ALLOCATED(NodeXValue)) ALLOCATE(NodeXValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeYValue)) ALLOCATE(NodeYValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeZValue)) ALLOCATE(NodeZValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeUValue)) ALLOCATE(NodeUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValue)) ALLOCATE(NodeVValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValue)) ALLOCATE(NodeWValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeUValueORG)) ALLOCATE(NodeUValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeVValueORG)) ALLOCATE(NodeVValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeWValueORG)) ALLOCATE(NodeWValueORG(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePValue)) ALLOCATE(NodePValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeMUValue)) ALLOCATE(NodeMUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeLabelValue)) ALLOCATE(NodeLabelValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeRHOValue)) ALLOCATE(NodeRHOValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeKappaValue)) ALLOCATE(NodeKappaValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,NodesPerElement(1)))
    IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,NodesPerElement(1)))

    CALL ENTERS("CMGUI OUTPUT",ERR,ERROR,*999)

    FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field

    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATED_POINT)
    
    CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)

    DO I=1,NumberOfElements

      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters

!       DO J=1,NodesPerElement(1)
      DO J=1,NodesPerElement(I)

        ELEMENT_NUMBER=I
! ! !         XI_COORDINATES(1)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
! ! !           & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,1)-1.0)/(REGION% &
! ! !           & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
! ! !           & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xi(1)-1.0)
! ! !         XI_COORDINATES(2)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
! ! !           & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,2)-1.0)/(REGION% &
! ! !           & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
! ! !           & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xi(2)-1.0)
! ! !         XI_COORDINATES(3)=(REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
! ! !           & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%node_position_index(J,3)-1.0)/(REGION% &
! ! !           & equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations%interpolation% &
! ! !           & geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%number_of_nodes_xi(3)-1.0)
! ! !         IF(NumberOfDimensions==2)THEN
! ! !           STOP 'Encas format only available for 3D hex and tets'
! ! !         END IF
! ! ! 
        !K is global node number
        K=REGION%meshes%meshes(1)%ptr%topology(1)%ptr%elements%elements(I)%global_element_nodes(J)

        COORDINATES=(/1,1,1/)
! ! ! 
! ! !         CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
! ! !           & INTERPOLATION_PARAMETERS(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
! ! !         CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,INTERPOLATED_POINT(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
        NodeXValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)
        NodeYValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))

        IF(NumberOfDimensions==3)THEN
          NodeZValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
        END IF

        IF((EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_POISSON_EQUATION_TYPE) &
          & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE))THEN 
          NodeUValue(K)=0.0_DP
          NodeVValue(K)=0.0_DP
        ELSE
!         NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field%variables(1) &
          NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field% &
            & variables(var_idx)%parameter_sets%parameter_sets(2)%ptr%parameters%cmiss%data_dp(K)
! !         NodeVValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field%variables(1) &
          NodeVValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field% &
            & variables(var_idx)%parameter_sets%parameter_sets(2)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))
        ENDIF

        parameter_set_idx = 2
        NodeUValueORG(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
          & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)
        parameter_set_idx = 3
        NodeVValueORG(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
          & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)

        IF(NumberOfDimensions==3)THEN
          parameter_set_idx = 4
          NodeWValueORG(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
            & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)

          IF((EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_POISSON_EQUATION_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE))THEN 
            NodeWValue(K)=0.0_DP
          ELSE! 
! !           NodeWValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field%variables(1) &
            NodeWValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field% &
              & variables(var_idx)%parameter_sets%parameter_sets(2)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
          END IF
        END IF


! ! !       NodeUValue(K)=INTERPOLATED_POINT%VALUES(1,1)
! ! !       NodeVValue(K)=INTERPOLATED_POINT%VALUES(2,1)
! ! !       NodeWValue(K)=INTERPOLATED_POINT%VALUES(3,1)

        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) &
          & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)   &
          & .OR. (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_POISSON_EQUATION_TYPE) &
              & .AND.((EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE) &
              & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) &
                & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE))) THEN
          parameter_set_idx = 1 
          NodePValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
!             & variables(1)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)
            & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)
        END IF

          NodeMUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1)
          NodeRHOValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(2)

          parameter_set_idx = 5
          NodeLabelValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
            & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)

      END DO 
    END DO

!     NodeMUValue=0.0_DP
!     NodeRHOValue=0.0_DP

    IF( NumberOfDimensions==3 )THEN
      !For 3D, the following call works ...
      lagrange_simplex=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations% &
        & interpolation%geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%type
!         lagrange_simplex=2
    ELSE
      !chrm, 20.08.09:
      ! ... but the above call does not work for 2D.
      !Thus, for 2D, we hard-wire it to 'quad':
      IF(MaxNodesPerElement==4.OR.MaxNodesPerElement==9.OR.MaxNodesPerElement==16) THEN
        lagrange_simplex=1
      ELSE IF(MaxNodesPerElement==3.OR.MaxNodesPerElement==6.OR.MaxNodesPerElement==10) THEN
        lagrange_simplex=2
      ENDIF
    END IF

    !This is for Poisson-Flow problems only
    IF(NumberOfVariableComponents==1)NumberOfVariableComponents=NumberOfDimensions+1
    IF(NumberOfMaterialComponents==5)NumberOfMaterialComponents=3
    IF(NumberOfMaterialComponents==4)NumberOfMaterialComponents=3
    IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS &
      & .AND.EQUATIONS_SET%TYPE==EQUATIONS_SET_POISSON_EQUATION_TYPE &
      & .AND.EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
         NumberOfVariableComponents=NumberOfVariableComponents+NumberOfDimensions
    ENDIF 


    NumberOfFieldComponent(1)=NumberOfDimensions
    NumberOfFieldComponent(2)=NumberOfVariableComponents
    NumberOfFieldComponent(3)=NumberOfMaterialComponents

    DO I=1,NumberOfElements
!       DO J=1,NodesPerElement(1)
      DO J=1,NodesPerElement(I)
        ElementNodes(I,J)=REGION%meshes%meshes(1)%ptr%topology(1)% &
          & ptr%elements%elements(I)%global_element_nodes(J)
        ElementNodesScales(I,J)=1.0000000000000000E+00
      END DO
    END DO

    CALL FLUID_MECHANICS_IO_WRITE_NODES_CMGUI(NAME,EQUATIONS_SET)
    CALL FLUID_MECHANICS_IO_WRITE_ELEMENTS_CMGUI(NAME)
    CALL FLUID_MECHANICS_IO_WRITE_DATA_ENCAS_BLOCK(NAME)


    CALL EXITS("FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK",ERR,ERROR)    
    CALL EXITS("FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK")
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK


  ! OK
  !================================================================================================================================
  !  

  !> Executes nodes writing process.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_DATA_ENCAS_BLOCK(NAME)

    IMPLICIT NONE
    CHARACTER(14), INTENT(IN) :: NAME !<the prefix name of file.
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
!     CHARACTER :: FILENAME !<the prefix name of file.
    INTEGER:: I
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR 
!     DOUBLE PRECISION:: velocity_magnitude
!     CHARACTER(80):: OUTPUT_FILE

  !==================

    FILENAME="./output/"//NAME//".scl1"
    OPEN(UNIT=81, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(81,*)'Absolute Pressure' 
    WRITE(81,*)'part'
    WRITE(81,*)'        1'
    WRITE(81,*)'block'
    DO I = 1,NodesPerMeshComponent(1)
      WRITE(81,'(e12.5)')NodePValue(I)
    ENDDO
    CLOSE(81)
 
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Encas data...",ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_DATA_ENCAS_BLOCK",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_DATA_ENCAS_BLOCK



  ! OK
  !================================================================================================================================
  !

  !> Executes nodes writing process.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_NODES_CMGUI(NAME,EQUATIONS_SET)

    IMPLICIT NONE

    CHARACTER(14), INTENT(IN) :: NAME !<the prefix name of file.
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
!     CHARACTER :: FILENAME !<the prefix name of file.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG):: I
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR 
    LOGICAL:: ANALYTIC

    ANALYTIC=.FALSE.
    IF( DARCY%ANALYTIC ) ANALYTIC=.TRUE.
 
    IF( DARCY%ANALYTIC ) THEN
      CALL FLUID_MECHANICS_IO_DARCY_GET_ANALYTIC
      CALL FLUID_MECHANICS_IO_DARCY_EVAL_ERROR
    END IF
 
    FILENAME="./output/"//NAME//".exnode"
    OPEN(UNIT=14, FILE=CHAR(FILENAME),STATUS='unknown')

! WRITING HEADER INFORMATION

    WRITE(14,*) 'Group name: OpenCMISS'

    IF( ANALYTIC ) THEN
      WRITE(14,*) '#Fields=',TRIM(NMs(NumberOfFields + 2))
    ELSE
      WRITE(14,*) '#Fields=',TRIM(NMs(NumberOfFields))
    END IF

    ValueIndex=1
    WRITE(14,*) ' 1) coordinates,  coordinate, rectangular cartesian, #Components=',TRIM(NMs(NumberOfDimensions))

    DO I=1,NumberOfDimensions
      IF(I==1) THEN
        WRITE(14,*) '   x.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
      ELSE IF(I==2) THEN
        WRITE(14,*) '   y.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
      ELSE
        WRITE(14,*) '   z.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
      END IF
      ValueIndex=ValueIndex+1
    END DO

    WRITE(14,*) ' 2) general,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))

    DO I=1,NumberOfVariableComponents
      WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
      ValueIndex=ValueIndex+1
    END DO

    WRITE(14,*) ' 3) material,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfMaterialComponents))

    DO I=1,NumberOfMaterialComponents
      WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
      ValueIndex=ValueIndex+1
    END DO

    IF( ANALYTIC ) THEN
      WRITE(14,*) ' 4) exact,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))
      DO I=1,NumberOfVariableComponents
        WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
        ValueIndex=ValueIndex+1
      END DO

      WRITE(14,*) ' 5) error,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))
      DO I=1,NumberOfVariableComponents
        WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
        ValueIndex=ValueIndex+1
      END DO
    END IF


! NOW WRITE NODE INFORMATION

    DO I = 1,NodesPerMeshComponent(1)
      WRITE(14,*) ' Node: ',I
      WRITE(14,'("    ", es25.16 )')NodeXValue(I)

      IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
        WRITE(14,'("    ", es25.16 )')NodeYValue(I)
      END IF

      IF(NumberOfDimensions==3) THEN
        WRITE(14,'("    ", es25.16 )')NodeZValue(I)
      END IF

      WRITE(14,'("    ", es25.16 )')NodeUValue(I)
      WRITE(14,'("    ", es25.16 )')NodeVValue(I)

      IF(NumberOfDimensions==3 .OR. NumberOfDimensions==1) THEN
      WRITE(14,'("    ", es25.16 )')NodeWValue(I)
      END IF

      IF(NumberOfDimensions/=1) THEN
        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) &
          & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
                & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
                & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE))THEN
            WRITE(14,'("    ", es25.16 )')NodePValue(I)
        END IF
        IF((EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_POISSON_EQUATION_TYPE) &
          & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)) THEN
          WRITE(14,'("    ", es25.16 )')NodeUValueORG(I)
          WRITE(14,'("    ", es25.16 )')NodeVValueORG(I)
          WRITE(14,'("    ", es25.16 )')NodeWValueORG(I)
        END IF
        IF((EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_POISSON_EQUATION_TYPE) &
          & .AND.((EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE) &
          & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) &
          & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE))) THEN
          WRITE(14,'("    ", es25.16 )')NodePValue(I)
          WRITE(14,'("    ", es25.16 )')NodeLabelValue(I)
        END IF
      END IF

!---tob: Mass increase for coupled elasticity Darcy INRIA model
        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
          & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE) &
            & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
              IF(NumberOfDimensions==3)THEN
                WRITE(14,'("    ", es25.16 )')NodeMIValue(I)
              END IF
        END IF
!---toe

      WRITE(14,'("    ", es25.16 )')NodeMUValue(I)
      WRITE(14,'("    ", es25.16 )')NodeRHOValue(I)

      IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
        & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_1DTRANSIENT_NAVIER_STOKES_SUBTYPE)) THEN
        WRITE(14,'("    ", es25.16 )')NodeEValue(I)
        WRITE(14,'("    ", es25.16 )')NodeH0Value(I)
        WRITE(14,'("    ", es25.16 )')NodeA0Value(I)
        WRITE(14,'("    ", es25.16 )')NodeSIGMAValue(I)
      END IF

      IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS)THEN
        IF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE)THEN
          IF( (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
            & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) &
            & .OR.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) &
            & .OR. (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) )THEN
            WRITE(14,'("    ", es25.16 )')NodeKappaValue(I)
          END IF
        END IF
      END IF

      IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) & 
        & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_DARCY_EQUATION_TYPE) &
          & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) )THEN
          !--- Remaining tensor material data of permeability tensor
          WRITE(14,'("    ", es25.16 )')NodePerm2Value(I)
          WRITE(14,'("    ", es25.16 )')NodePerm3Value(I)
          WRITE(14,'("    ", es25.16 )')NodePerm4Value(I)
          WRITE(14,'("    ", es25.16 )')NodePerm5Value(I)
          WRITE(14,'("    ", es25.16 )')NodePerm6Value(I)
      END IF


      IF( ANALYTIC ) THEN
        WRITE(14,'("    ", es25.16 )')NodeUValue_analytic(I)
        IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
          WRITE(14,'("    ", es25.16 )')NodeVValue_analytic(I)
        END IF
        IF(NumberOfDimensions==3) THEN
          WRITE(14,'("    ", es25.16 )')NodeWValue_analytic(I)
        END IF
        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) &
          & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
          WRITE(14,'("    ", es25.16 )')NodePValue_analytic(I)
        END IF

        WRITE(14,'("    ", es25.16 )')NodeUValue_error(I)
        IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
          WRITE(14,'("    ", es25.16 )')NodeVValue_error(I)
        END IF
        IF(NumberOfDimensions==3) THEN
          WRITE(14,'("    ", es25.16 )')NodeWValue_error(I)
        END IF
        IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS) &
          & .OR.(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS) & 
              & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE.NE.EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) )THEN
          WRITE(14,'("    ", es25.16 )')NodePValue_error(I)
        END IF
      END IF

    END DO
 
    WRITE(14,*) ' '
    CLOSE(14)

    IF( DARCY%ANALYTIC ) THEN
      CALL FLUID_MECHANICS_IO_DARCY_EVAL_MAX_ERROR
    END IF

!test output for DN only
    IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
      IF(DN) THEN
        FILENAME="./output/"//NAME//".davidn"
        OPEN(UNIT=14, FILE=CHAR(FILENAME),STATUS='unknown')
        WRITE(14,*) NodesPerMeshComponent(1),NodesPerMeshComponent(1),NodesPerMeshComponent(2)
        DO I=1,NodesPerMeshComponent(1) 
          WRITE(14,'(3("    ", es25.16 ))')NodeXValue(I),NodeYValue(I),NodeZValue(I)
        ENDDO
        DO I=1,NodesPerMeshComponent(1) 
          WRITE(14,'(6("    ", es25.16 ))')NodeXValue(I),NodeYValue(I),NodeZValue(I),NodeUValue(I),NodeVValue(I),NodeWValue(I)
        ENDDO
        DO I=1,NodesPerMeshComponent(2)
          WRITE(14,'(6("    ", es25.16 ))')NodeXValue(I),NodeYValue(I),NodeZValue(I),NodePValue2(I)
        ENDDO
        CLOSE(14)
      ENDIF
    END IF

    IF(NumberOfDimensions==1) THEN
      IF(DN) THEN
        FILENAME="./output/"//NAME//".davidn"
        OPEN(UNIT=14, FILE=CHAR(FILENAME),STATUS='unknown')
        WRITE(14,*) NodesPerMeshComponent(1),NodesPerMeshComponent(1)
        DO I=1,NodesPerMeshComponent(1) 
          WRITE(14,'(3("    ", es25.16 ))')NodeXValue(I),NodeYValue(I),NodeZValue(I)
        ENDDO
        DO I=1,NodesPerMeshComponent(1) 
          WRITE(14,'(6("    ", es25.16 ))')NodeXValue(I),NodeYValue(I),NodeZValue(I),NodeUValue(I),NodeVValue(I),NodeWValue(I)
        ENDDO
        CLOSE(14)
      ENDIF
    END IF

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Nodes...",ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_NODES_CMGUI",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_NODES_CMGUI

  ! OK
  !================================================================================================================================
  !

  !> Executes element writing process.
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_ELEMENTS_CMGUI(NAME)

!     TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    CHARACTER(14), INTENT(IN) :: NAME !<the prefix name of file.
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
!     CHARACTER :: FILENAME !<the prefix name of file.
    ! CHARACTER*60 ELEM_TYPE
    INTEGER(INTG):: I,J,K,KK
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR
    LOGICAL:: OUTPUT_FLAG

    FILENAME="./output/"//NAME//".exelem"
    OPEN(UNIT=5, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(5,*) 'Group name: OpenCMISS'


  DO KK = 1,NumberOfElements

    !write out the running header only once #NodesPerElement is changing
    IF( KK == 1 ) THEN
      OUTPUT_FLAG = .TRUE.
    ELSE IF( KK > 1 .AND. NodesPerElement(KK).NE.NodesPerElement(KK-1) ) THEN
      OUTPUT_FLAG = .TRUE.
    ELSE 
      OUTPUT_FLAG = .FALSE.
    END IF

    IF( OUTPUT_FLAG ) THEN
 
    IF(lagrange_simplex==2) THEN
      IF(NumberOfDimensions==2) THEN
        WRITE(5,*) 'Shape.  Dimension=',TRIM(NMs(NumberOfDimensions)),', simplex(2)*simplex'
        IF(MaxNodesPerElement==3) THEN
          WRITE(5,*) '#Scale factor sets= 1'
!           WRITE(5,*) ' l.simplex(2)*l.simplex, #Scale factors= ', NodesPerElement(1)
          WRITE(5,*) ' l.simplex(2)*l.simplex, #Scale factors= ', NodesPerElement(KK)
        ELSE IF(MaxNodesPerElement==6) THEN
          WRITE(5,*) '#Scale factor sets= 1'
!           WRITE(5,*) ' l.simplex(2)*l.simplex, #Scale factors= ', NodesPerElement(1)
          WRITE(5,*) ' l.simplex(2)*l.simplex, #Scale factors= ', NodesPerElement(KK)
        ELSE IF (MaxNodesPerElement== 10 ) THEN
          WRITE(5,*) '#Scale factor sets= 1'
!           WRITE(5,*) ' q.simplex(2)*q.simplex, #Scale factors= ', NodesPerElement(1)
          WRITE(5,*) ' q.simplex(2)*q.simplex, #Scale factors= ', NodesPerElement(KK)
        ENDIF
      ELSE IF(NumberOfDimensions==3) THEN
        WRITE(5,*) 'Shape.  Dimension=',TRIM(NMs(NumberOfDimensions)),', simplex(2;3)*simplex*simplex'
        IF(MaxNodesPerElement==4) THEN
          WRITE(5,*) '#Scale factor sets= 1'
!           WRITE(5,*) ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', NodesPerElement(1)
          WRITE(5,*) ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', NodesPerElement(KK)
        ELSE IF (MaxNodesPerElement== 10 ) THEN
          WRITE(5,*) '#Scale factor sets= 1'
!           WRITE(5,*) ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', NodesPerElement(1)
          WRITE(5,*) ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', NodesPerElement(KK)
        ELSE IF(MaxNodesPerElement==20) THEN
          WRITE(5,*) '#Scale factor sets= 1'
!           WRITE(5,*) ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', NodesPerElement(1)
          WRITE(5,*) ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', NodesPerElement(KK)
        ENDIF      
      ELSE
        WRITE(5,*) '#Scale factor sets= 0'
      END IF
    ELSE IF (lagrange_simplex==1) THEN
      WRITE(5,*) 'Shape.  Dimension= ',TRIM(NMs(NumberOfDimensions))
      WRITE(5,*) '#Scale factor sets= 1'
      IF(NumberOfDimensions==1) THEN
!            WRITE(5,*) 'q.Lagrange, #Scale factors=',NodesPerElement(1)
           WRITE(5,*) 'q.Lagrange, #Scale factors=',NodesPerElement(KK)
      ELSE IF (NumberOfDimensions==2) THEN
        IF(MaxNodesPerElement==4) THEN
!           WRITE(5,*) 'l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(1)
          WRITE(5,*) 'l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(KK)
        ELSE IF(MaxNodesPerElement==9) THEN
!           WRITE(5,*) 'q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(1)
          WRITE(5,*) 'q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(KK)
        ELSE IF(MaxNodesPerElement==16) THEN
!           WRITE(5,*) 'c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(1)
          WRITE(5,*) 'c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(KK)
        END IF
      ELSE
        IF(MaxNodesPerElement==8) THEN
!           WRITE(5,*) 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(1)
          WRITE(5,*) 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(KK)
        ELSE IF(MaxNodesPerElement==27) THEN
!           WRITE(5,*) 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(1)
          WRITE(5,*) 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(KK)
        ELSE IF(MaxNodesPerElement==64) THEN
!           WRITE(5,*) 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(1)
          WRITE(5,*) 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(KK)
        END IF
      END IF
    END IF

!     WRITE(5,*) '#Nodes= ',TRIM(NMs(NodesPerElement(1)))
    WRITE(5,*) '#Nodes= ',TRIM(NMs(NodesPerElement(KK)))
    WRITE(5,*) '#Fields= ',TRIM(Nms(NumberOfFields))

    DO I=1,NumberOfFields
      IF(I==1)THEN
        WRITE(5,*)' 1) coordinates,  coordinate, rectangular cartesian, #Components= ',TRIM(NMs(NumberOfDimensions))
      ELSE IF(I==2) THEN
        WRITE(5,*)' 2) general,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfVariableComponents))
      ELSE IF(I==3) THEN
        WRITE(5,*)' 3) material,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfMaterialComponents))
      END IF

      DO J=1,NumberOfFieldComponent(I)
        IF(NumberOfDimensions==1) THEN
          IF(I==1)THEN
            IF(J==1) THEN
                WRITE(5,*)'   x.   q.Lagrange, no modify, standard node based.'
            ELSE IF(J==2) THEN
                WRITE(5,*)'   y.   q.Lagrange, no modify, standard node based.'
            ELSE IF(J==3) THEN
                WRITE(5,*)'   z.   q.Lagrange, no modify, standard node based.'
            END IF
          ELSE
            WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange, no modify, standard node based.'
          END IF
        ELSE IF(NumberOfDimensions==2) THEN
          IF(I==1)THEN
            IF(J==1) THEN
              IF(MaxNodesPerElement==4)THEN
                WRITE(5,*)'   x.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9) THEN
                WRITE(5,*)'   x.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(5,*)'   x.   c.Lagrange*c.Lagrange, no modify, standard node based.'

              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(5,*)'   x.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(5,*)'   x.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   x.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF 
            ELSE IF(J==2) THEN
              IF(MaxNodesPerElement==4) THEN
                WRITE(5,*)'   y.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(5,*)'   y.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(5,*)'   y.   c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(5,*)'   y.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(5,*)'   y.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   y.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF
            ELSE IF(J==3) THEN
              IF(MaxNodesPerElement==4) THEN
                WRITE(5,*)'   z.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(5,*)'   z.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(5,*)'   z.   c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(5,*)'   z.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(5,*)'   z.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   z.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF
            END IF
          ELSE
              IF(MaxNodesPerElement==4) THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF
          END IF
        ELSE IF(NumberOfDimensions==3) THEN
          IF(I==1)THEN
            IF(J==1) THEN
              IF(MaxNodesPerElement==8) THEN
                WRITE(5,*)'   x.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(5,*)'   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(5,*)'   x.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(5,*)'   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(5,*)'   x.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF 
            ELSE IF(J==2) THEN
              IF(MaxNodesPerElement==8) THEN
                WRITE(5,*)'   y.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(5,*)'   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(5,*)'   y.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(5,*)'   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   y.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(5,*)'   y.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF
            ELSE IF(J==3) THEN
              IF(MaxNodesPerElement==8) THEN
                WRITE(5,*)'   z.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(5,*)'   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(5,*)'   z.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(5,*)'   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   z.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(5,*)'   z.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF
            END IF
          ELSE
              IF(MaxNodesPerElement==8) THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF
          END IF
        END IF

        WRITE(5,*) '   #Nodes= ',TRIM(NMs(MaxNodesPerElement))
 
        DO K = 1,MaxNodesPerElement
          WRITE(5,*) '    ',TRIM(NMs(K)),'.  #Values=1'
          WRITE(5,*) '     Value indices:     1'
          WRITE(5,*) '     Scale factor indices:   ',TRIM(NMs(K))
        END DO
      END DO
    END DO

    END IF !write out the running header only once #NodesPerElement is changing


    IF(lagrange_simplex==2) THEN
!       IF(.NOT.ALLOCATED(SimplexOutputHelp)) ALLOCATE(SimplexOutputHelp(NodesPerElement(1)))
      IF(.NOT.ALLOCATED(SimplexOutputHelp)) ALLOCATE(SimplexOutputHelp(NodesPerElement(KK)))

!       DO K = 1,NumberOfElements
      K = KK

        IF(NumberOfDimensions==2)THEN
          SimplexOutputHelp(1)=ElementNodes(K,1)
          SimplexOutputHelp(2)=ElementNodes(K,4)
          SimplexOutputHelp(3)=ElementNodes(K,2)
          SimplexOutputHelp(4)=ElementNodes(K,6)
          SimplexOutputHelp(5)=ElementNodes(K,5)
          SimplexOutputHelp(6)=ElementNodes(K,3)
        ELSE IF(NumberOfDimensions==3) THEN
          SimplexOutputHelp(1)=ElementNodes(K,1)
          SimplexOutputHelp(2)=ElementNodes(K,5)
          SimplexOutputHelp(3)=ElementNodes(K,2)
          SimplexOutputHelp(4)=ElementNodes(K,7)
          SimplexOutputHelp(5)=ElementNodes(K,10)
          SimplexOutputHelp(6)=ElementNodes(K,4)
          SimplexOutputHelp(7)=ElementNodes(K,6)
          SimplexOutputHelp(8)=ElementNodes(K,8)
          SimplexOutputHelp(9)=ElementNodes(K,9)
          SimplexOutputHelp(10)=ElementNodes(K,3)
        END IF

        WRITE(5,*) 'Element:     ', K,' 0  0'
        WRITE(5,*) '   Nodes:'
        WRITE(5,*) '   ', SimplexOutputHelp
        WRITE(5,*) '   Scale factors:'
!         WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(1))
        WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(K))
!       END DO

    ELSE IF (lagrange_simplex==1) THEN

!       DO K = 1,NumberOfElements
      K = KK

        WRITE(5,*) 'Element:     ', K,' 0  0'
        WRITE(5,*) '   Nodes:'
!         WRITE(5,*) '   ', ElementNodes(K,1:NodesPerElement(1))
        WRITE(5,*) '   ', ElementNodes(K,1:NodesPerElement(K))
        WRITE(5,*) '   Scale factors:'
!         WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(1))
        WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(K))
!       END DO
    END IF


  ENDDO



    WRITE(5,*) ' '
    CLOSE(5)
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Elements...",ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_WRITE_ELEMENTS_CMGUI",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_ELEMENTS_CMGUI

  ! OK
  !================================================================================================================================
  !

  !> Reads in information defined by cmheart input file format.
  SUBROUTINE FLUID_MECHANICS_IO_READ_CMHEART1(EXPORT,ERR)

    !Argument variables
    TYPE (EXPORT_CONTAINER):: EXPORT  
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING):: ERROR !<The error string
    !Local Variables
    TYPE (EXPORT_CONTAINER):: TMP  
    ! INTEGER(INTG):: test 

    CALL ENTERS("FLUID_MECHANICS_IO_READ_CMHEART1",ERR,ERROR,*999)
    WRITE(*,*)' '
    WRITE(*,*)'Importing CMHEART information...'
    WRITE(*,*)' '
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Run from bin directory if input needed.",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"REMEMBER: M,V,P need to be defined as required by cmHeart!",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
!     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Press ENTER to start.",ERR,ERROR,*999)
!     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
!     READ(*,*)
    OPEN(UNIT=42, FILE='./input/CMHEART.inp',STATUS='old')

    CALL FLUID_MECHANICS_IO_READ_AUX
    CALL FLUID_MECHANICS_IO_READ_NODES
    CALL FLUID_MECHANICS_IO_READ_ELEMENTS

    IF(.NOT.ALLOCATED(OPENCMISS_ELEM_M)) ALLOCATE(OPENCMISS_ELEM_M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)), & 
      & STAT=ALLOC_ERROR)
    IF(.NOT.ALLOCATED(OPENCMISS_ELEM_V))ALLOCATE(OPENCMISS_ELEM_V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)), &
      & STAT=ALLOC_ERROR)
    IF(.NOT.ALLOCATED(OPENCMISS_ELEM_P))ALLOCATE(OPENCMISS_ELEM_P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)), &
      & STAT=ALLOC_ERROR)

    CALL FLUID_MECHANICS_IO_MAKE_UNIQUE
    CALL FLUID_MECHANICS_IO_ORDER_NUMBERING(OPENCMISS_ELEM_M,MESH_INFO(1)%T,NumberOfElementsDefined(1), & 
      & NumberOfNodesPerElement(1),1)
    CALL FLUID_MECHANICS_IO_ORDER_NUMBERING(OPENCMISS_ELEM_V,MESH_INFO(2)%T,NumberOfElementsDefined(2), & 
      & NumberOfNodesPerElement(2),2)
    CALL FLUID_MECHANICS_IO_ORDER_NUMBERING(OPENCMISS_ELEM_P,MESH_INFO(3)%T,NumberOfElementsDefined(3), & 
      & NumberOfNodesPerElement(3),3)
    WRITE(*,*)' '
    WRITE(*,*)'Import finished successfully...'
    WRITE(*,*)' '
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Export finished successfully...",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)

    IF(ALLOC_ERROR.NE.0) THEN
      CALL FLAG_ERROR("Error during allocation.",ERR,ERROR,*999)
    END IF

    ALLOCATE(TMP%M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)),STAT=ALLOC_ERROR)
    ALLOCATE(TMP%V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)),STAT=ALLOC_ERROR)
    ALLOCATE(TMP%P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)),STAT=ALLOC_ERROR)
    ALLOCATE(TMP%N(TotalNumberOfNodes,3),STAT=ALLOC_ERROR)

    TMP%M=OPENCMISS_ELEM_M
    TMP%V=OPENCMISS_ELEM_V
    TMP%P=OPENCMISS_ELEM_P
    TMP%N=OPENCMISS_NODE_COORD
    TMP%D=DIMEN
    TMP%F=BASE_INFO%n_B
    TMP%ID_M=1
    TMP%ID_V=2
    TMP%ID_P=3
    TMP%IT_M=OPENCMISS_INTERPOLATION(1)
    TMP%IT_V=OPENCMISS_INTERPOLATION(2)
    TMP%IT_P=OPENCMISS_INTERPOLATION(3)
  
    IF (BASE_INFO%HEXA==1) THEN
    !LAGRANGIAN BASIS
      TMP%IT_T=1
    ELSE 
    ! SIMPLEX BASIS
      TMP%IT_T=2
    END IF

    TMP%E_M=NumberOfElementsDefined(1)
    TMP%E_V=NumberOfElementsDefined(2)
    TMP%E_P=NumberOfElementsDefined(3)
    TMP%E_T=NumberOfElementsDefined(3)
    TMP%EN_M=NumberOfNodesPerElement(1)
    TMP%EN_V=NumberOfNodesPerElement(2)
    TMP%EN_P=NumberOfNodesPerElement(3)
    TMP%EN_T=TMP%EN_M+TMP%EN_V+TMP%EN_P
    TMP%N_M=ArrayOfNodesDefined(1)
    TMP%N_V=ArrayOfNodesDefined(2)
    TMP%N_P=ArrayOfNodesDefined(3)
    TMP%N_T=TMP%N_M+TMP%N_V+TMP%N_P

    EXPORT=TMP

    IF(ALLOC_ERROR.NE.0) THEN
      CALL FLAG_ERROR("Error during allocation.",ERR,ERROR,*999)
    END IF

    CALL EXITS("FLUID_MECHANICS_IO_READ_CMHEART1")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_READ_CMHEART1",ERR,ERROR)    
    CALL EXITS("FLUID_MECHANICS_IO_READ_CMHEART1")
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_CMHEART1


  ! OK
  !================================================================================================================================
  !


  !> Reads in information defined by cmheart input file format.
  SUBROUTINE FLUID_MECHANICS_IO_READ_CMHEART3(EXPORT1,EXPORT2,EXPORT3,CONNECT,ERR)

    !Argument variables
    TYPE (EXPORT_CONTAINER):: EXPORT1,EXPORT2,EXPORT3
    TYPE (COUPLING_PARAMETERS):: CONNECT
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING):: ERROR !<The error string
    !Local Variables
    TYPE (EXPORT_CONTAINER):: TMP1  
    INTEGER(INTG):: I

    CALL ENTERS("FLUID_MECHANICS_IO_READ_CMHEART3",ERR,ERROR,*999)
    WRITE(*,*)' '
    WRITE(*,*)'Importing CMHEART information...'
    WRITE(*,*)' '
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Run from bin directory if input needed.",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"REMEMBER: M,V,P need to be defined as required by cmHeart!",ERR,ERROR,*999)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
!     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Press ENTER to start.",ERR,ERROR,*999)
!     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)
!     READ(*,*)

    DO I=1,3
      IF(I==1) THEN
        OPEN(UNIT=42, FILE='./input/CMHEART1.inp',STATUS='old')
      ELSE IF(I==2) THEN
        OPEN(UNIT=42, FILE='./input/CMHEART2.inp',STATUS='old')
      ELSE IF(I==3) THEN
        OPEN(UNIT=42, FILE='./input/CMHEART3.inp',STATUS='old')
      ENDIF

      CALL FLUID_MECHANICS_IO_READ_AUX
      CALL FLUID_MECHANICS_IO_READ_NODES
      CALL FLUID_MECHANICS_IO_READ_ELEMENTS

      IF(ALLOCATED(OPENCMISS_ELEM_M)) DEALLOCATE(OPENCMISS_ELEM_M)
      IF(ALLOCATED(OPENCMISS_ELEM_V)) DEALLOCATE(OPENCMISS_ELEM_V)
      IF(ALLOCATED(OPENCMISS_ELEM_P)) DEALLOCATE(OPENCMISS_ELEM_P)

      IF(.NOT.ALLOCATED(OPENCMISS_ELEM_M)) ALLOCATE(OPENCMISS_ELEM_M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)), & 
        & STAT=ALLOC_ERROR)
      IF(.NOT.ALLOCATED(OPENCMISS_ELEM_V))ALLOCATE(OPENCMISS_ELEM_V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)), &
        & STAT=ALLOC_ERROR)
      IF(.NOT.ALLOCATED(OPENCMISS_ELEM_P))ALLOCATE(OPENCMISS_ELEM_P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)), &
        & STAT=ALLOC_ERROR)

      CALL FLUID_MECHANICS_IO_MAKE_UNIQUE
      CALL FLUID_MECHANICS_IO_ORDER_NUMBERING(OPENCMISS_ELEM_M,MESH_INFO(1)%T,NumberOfElementsDefined(1), & 
        & NumberOfNodesPerElement(1),1)
      CALL FLUID_MECHANICS_IO_ORDER_NUMBERING(OPENCMISS_ELEM_V,MESH_INFO(2)%T,NumberOfElementsDefined(2), & 
        & NumberOfNodesPerElement(2),2)
      CALL FLUID_MECHANICS_IO_ORDER_NUMBERING(OPENCMISS_ELEM_P,MESH_INFO(3)%T,NumberOfElementsDefined(3), & 
        & NumberOfNodesPerElement(3),3)
      WRITE(*,*)' '
      WRITE(*,*)'Import finished successfully...', I
      WRITE(*,*)' '
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Export finished successfully...",ERR,ERROR,*999)
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," ",ERR,ERROR,*999)

      IF(ALLOC_ERROR.NE.0) THEN
        CALL FLAG_ERROR("Error during allocation.",ERR,ERROR,*999)
      END IF

      ALLOCATE(TMP1%M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)),STAT=ALLOC_ERROR)
      ALLOCATE(TMP1%V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)),STAT=ALLOC_ERROR)
      ALLOCATE(TMP1%P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)),STAT=ALLOC_ERROR)
      ALLOCATE(TMP1%N(TotalNumberOfNodes,3),STAT=ALLOC_ERROR)

      TMP1%M=OPENCMISS_ELEM_M
      TMP1%V=OPENCMISS_ELEM_V
      TMP1%P=OPENCMISS_ELEM_P
      TMP1%N=OPENCMISS_NODE_COORD
      TMP1%D=DIMEN
      TMP1%F=BASE_INFO%n_B
      TMP1%ID_M=1
      TMP1%ID_V=2
      TMP1%ID_P=3
      TMP1%IT_M=OPENCMISS_INTERPOLATION(1)
      TMP1%IT_V=OPENCMISS_INTERPOLATION(2)
      TMP1%IT_P=OPENCMISS_INTERPOLATION(3)
  
      IF (BASE_INFO%HEXA==1) THEN
      !LAGRANGIAN BASIS
        TMP1%IT_T=1
      ELSE 
      ! SIMPLEX BASIS
        TMP1%IT_T=2
      END IF

      TMP1%E_M=NumberOfElementsDefined(1)
      TMP1%E_V=NumberOfElementsDefined(2)
      TMP1%E_P=NumberOfElementsDefined(3)
      TMP1%E_T=NumberOfElementsDefined(3)
      TMP1%EN_M=NumberOfNodesPerElement(1)
      TMP1%EN_V=NumberOfNodesPerElement(2)
      TMP1%EN_P=NumberOfNodesPerElement(3)
      TMP1%EN_T=TMP1%EN_M+TMP1%EN_V+TMP1%EN_P
      TMP1%N_M=ArrayOfNodesDefined(1)
      TMP1%N_V=ArrayOfNodesDefined(2)
      TMP1%N_P=ArrayOfNodesDefined(3)
      TMP1%N_T=TMP1%N_M+TMP1%N_V+TMP1%N_P

      IF(I==1) THEN
        EXPORT1=TMP1
      ELSE IF (I==2) THEN
        EXPORT2=TMP1
      ELSE IF (I==3) THEN
        EXPORT3=TMP1
      ENDIF

      IF(ALLOC_ERROR.NE.0) THEN
        CALL FLAG_ERROR("Error during allocation.",ERR,ERROR,*999)
      END IF
    ENDDO


    !Read in the interface connectivity mapping
    OPEN(UNIT=79, FILE='./input/IM_COUPLING.LM',STATUS='old')
    READ(79,*) CONNECT%NUMBER_OF_COUPLINGS

    ALLOCATE(CONNECT%INTERFACE_ELEMENT_NUMBER(CONNECT%NUMBER_OF_COUPLINGS),STAT=ALLOC_ERROR)
    ALLOCATE(CONNECT%INTERFACE_ELEMENT_LOCAL_NODE(CONNECT%NUMBER_OF_COUPLINGS),STAT=ALLOC_ERROR)
    ALLOCATE(CONNECT%MESH1_ELEMENT_NUMBER(CONNECT%NUMBER_OF_COUPLINGS),STAT=ALLOC_ERROR)
    ALLOCATE(CONNECT%MESH1_ELEMENT_XI(CONNECT%NUMBER_OF_COUPLINGS,3),STAT=ALLOC_ERROR)
    ALLOCATE(CONNECT%MESH2_ELEMENT_NUMBER(CONNECT%NUMBER_OF_COUPLINGS),STAT=ALLOC_ERROR)
    ALLOCATE(CONNECT%MESH2_ELEMENT_XI(CONNECT%NUMBER_OF_COUPLINGS,3),STAT=ALLOC_ERROR)

    DO I=1,CONNECT%NUMBER_OF_COUPLINGS
      READ(79,*) CONNECT%INTERFACE_ELEMENT_NUMBER(I), &
        & CONNECT%INTERFACE_ELEMENT_LOCAL_NODE(I),CONNECT%MESH1_ID,CONNECT%MESH1_ELEMENT_NUMBER(I), &
        & CONNECT%MESH1_ELEMENT_XI(I,1:3), CONNECT%MESH2_ID,CONNECT%MESH2_ELEMENT_NUMBER(I), &
        & CONNECT%MESH2_ELEMENT_XI(I,1:3)
    ENDDO
    WRITE(*,*)' '
    WRITE(*,*)'Import finished successfully...  IMC'
    WRITE(*,*)' '

    IF(ALLOC_ERROR.NE.0) THEN
      CALL FLAG_ERROR("Error during allocation.",ERR,ERROR,*999)
    END IF

    CLOSE(79)

    CALL EXITS("FLUID_MECHANICS_IO_READ_CMHEART3")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_READ_CMHEART3",ERR,ERROR)    
    CALL EXITS("FLUID_MECHANICS_IO_READ_CMHEART3")
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_CMHEART3


  ! OK
  !================================================================================================================================
  !

  !> Executes cmheart parameter reading process.
  SUBROUTINE FLUID_MECHANICS_IO_READ_AUX

    INTEGER(INTG):: I
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR

    READ(42,*) NIMZ
    NIMZ = TRIM(NIMZ); BASE_INFO%n_B = 0; BASE_INFO%HEXA = 0; BASE_INFO%DM = 3
    OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read') ! Read base file for initial parameters

    BASE_INFO%n_B=0
    BASE_INFO%n_pts=0
    BASE_INFO%VL=0
    BASE_INFO%n_pts_f=0
    BASE_INFO%FACES=0
    BASE_INFO%FNODES=0
    BASE_INFO%HEXA=0
    BASE_INFO%DM=0
    BASE_INFO%TRI_BASIS=0
    BASE_INFO%TET_BASIS=0
    BASE_INFO%QUAD_BASIS=0
    BASE_INFO%HEX_BASIS=0


    DO WHILE (0 < 1)
      READ(1,*,END=50) IN_CHAR
      IF (INDEX(IN_CHAR,'no_fields!') == 1)         READ(1,*) BASE_INFO%n_B
      IF (INDEX(IN_CHAR,'no_gauss!') == 1)          READ(1,*) BASE_INFO%n_pts
      IF (INDEX(IN_CHAR,'volume!') == 1)            READ(1,*) BASE_INFO%VL
      IF (INDEX(IN_CHAR,'no_gauss_f!') == 1)        READ(1,*) BASE_INFO%n_pts_f
      IF (INDEX(IN_CHAR,'no_ele_faces!') == 1)      READ(1,*) BASE_INFO%FACES
      IF (INDEX(IN_CHAR,'no_ele_nodes_f!') == 1)    READ(1,*) BASE_INFO%FNODES
      IF (INDEX(IN_CHAR,'hexa_basis!') == 1)        BASE_INFO%HEXA = 1
      IF (INDEX(IN_CHAR,'domain_dimension!') == 1)  READ(1,*) BASE_INFO%DM
      IF (INDEX(IN_CHAR,'TRI_BASIS!') == 1)  BASE_INFO%TRI_BASIS = 1
      IF (INDEX(IN_CHAR,'TET_BASIS!') == 1)  BASE_INFO%TET_BASIS = 1
      IF (INDEX(IN_CHAR,'QUAD_BASIS!') == 1) BASE_INFO%QUAD_BASIS = 1
      IF (INDEX(IN_CHAR,'HEX_BASIS!') == 1)  BASE_INFO%HEX_BASIS = 1
    END DO

50  CLOSE(1)

    DIMEN=BASE_INFO%DM

    ALLOCATE(BASE_INFO%B(BASE_INFO%n_B),STAT=ALLOC_ERROR)
    OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read') ! Read base file for initial parameters

    DO WHILE (0 < 1)
      READ(1,*,END=52) IN_CHAR
      IF (INDEX(IN_CHAR,'no_basis_M!') == 1) READ(1,*) BASE_INFO%B(1)%n
      IF (INDEX(IN_CHAR,'no_basis_V!') == 1) READ(1,*) BASE_INFO%B(2)%n
      IF (INDEX(IN_CHAR,'no_basis_P!') == 1) READ(1,*) BASE_INFO%B(3)%n
      IF (INDEX(IN_CHAR,'dim_field_M!') == 1) READ(1,*) BASE_INFO%B(1)%DM
      IF (INDEX(IN_CHAR,'dim_field_V!') == 1) READ(1,*) BASE_INFO%B(2)%DM
      IF (INDEX(IN_CHAR,'dim_field_P!') == 1) READ(1,*) BASE_INFO%B(3)%DM
    END DO

52  CLOSE(1)

    IF (BASE_INFO%QUAD_BASIS/= 1.AND.BASE_INFO%HEX_BASIS /= 1.AND.BASE_INFO%TRI_BASIS/= 1.AND.BASE_INFO%TET_BASIS /= 1)THEN
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Cubic Hermite not implemented yet...",ERR,ERROR,*999)
    ELSE
 
      DO I=1,3
        IF (BASE_INFO%TRI_BASIS== 1.OR.BASE_INFO%TET_BASIS == 1)THEN
          IF(DIMEN==2) THEN
            IF(BASE_INFO%B(I)%n==3) THEN
              NumberOfNodesPerElement(I)=3
              OPENCMISS_INTERPOLATION(I)=7
            ELSE IF(BASE_INFO%B(I)%n==6) THEN
              NumberOfNodesPerElement(I)=6
              OPENCMISS_INTERPOLATION(I)=8
            ELSE IF(BASE_INFO%B(I)%n==10) THEN
              NumberOfNodesPerElement(I)=10
              OPENCMISS_INTERPOLATION(I)=9
            ELSE
              STOP
            END IF
          ELSE IF(DIMEN==3) THEN
            IF(BASE_INFO%B(I)%n==4) THEN
              NumberOfNodesPerElement(I)=4
              OPENCMISS_INTERPOLATION(I)=7
            ELSE IF(BASE_INFO%B(I)%n==10) THEN
              NumberOfNodesPerElement(I)=10
              OPENCMISS_INTERPOLATION(I)=8
            ELSE IF(BASE_INFO%B(I)%n==20) THEN
              NumberOfNodesPerElement(I)=20
              OPENCMISS_INTERPOLATION(I)=9
            ELSE
              STOP
            END IF
          ELSE 
            STOP
          END IF
        END IF
 
        IF (BASE_INFO%QUAD_BASIS== 1.OR.BASE_INFO%HEX_BASIS == 1)THEN
          IF(BASE_INFO%B(I)%n==2) THEN
            !2D/3D LINEAR LAGRANGE
            OPENCMISS_INTERPOLATION(I)=1
            IF(DIMEN==2) THEN
              NumberOfNodesPerElement(I)=4
            ELSE IF(DIMEN==3) THEN
              NumberOfNodesPerElement(I)=8
            ELSE 
              STOP
            END IF
          ELSE IF(BASE_INFO%B(I)%n==3) THEN
            !2D/3D QUADRATIC LAGRANGE
            OPENCMISS_INTERPOLATION(I)=2
            IF(DIMEN==2) THEN
              NumberOfNodesPerElement(I)=9
            ELSE IF(DIMEN==3) THEN
              NumberOfNodesPerElement(I)=27
            ELSE 
              STOP
            END IF
          ELSE IF(BASE_INFO%B(I)%n==4) THEN
            !2D/3D CUBIC LAGRANGE
            OPENCMISS_INTERPOLATION(I)=3
            IF(DIMEN==2) THEN
              NumberOfNodesPerElement(I)=16
            ELSE IF(DIMEN==3) THEN
              NumberOfNodesPerElement(I)=64
            ELSE 
              STOP
            END IF
          ELSE
            STOP
          END IF
        END IF
      END DO
    END IF

    IF(ALLOC_ERROR.NE.0) THEN
      STOP 'Error during allocation'
    END IF
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_READ_AUX",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_AUX

  ! OK
  !================================================================================================================================
  !

  !> Reorders the element node definition as needed by OpenCMISS.
  SUBROUTINE FLUID_MECHANICS_IO_ORDER_NUMBERING(NEW,OLD,n,m,I)

    INTEGER(INTG):: I,J,M,N
    INTEGER(INTG)::NEW(n,m),OLD(n,m)
!     INTEGER(INTG) :: ERR
!     TYPE(VARYING_STRING):: ERROR

!  NEW=OLD
    DO J=1,n
      IF (BASE_INFO%QUAD_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==2) THEN
          !2D HEX LINEAR
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==3) THEN
          !2D HEX QUADR
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,5)
          NEW(J,3)=OLD(J,2)
          NEW(J,4)=OLD(J,6)
          NEW(J,5)=OLD(J,7)
          NEW(J,6)=OLD(J,8)
          NEW(J,7)=OLD(J,3)
          NEW(J,8)=OLD(J,9)
          NEW(J,9)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==4) THEN
          !2D HEX CUB
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,5)
          NEW(J,3)=OLD(J,6)
          NEW(J,4)=OLD(J,2)
          NEW(J,5)=OLD(J,7)
          NEW(J,6)=OLD(J,8)
          NEW(J,7)=OLD(J,9)
          NEW(J,8)=OLD(J,10)
          NEW(J,9)=OLD(J,11)
          NEW(J,10)=OLD(J,12)
          NEW(J,11)=OLD(J,13)
          NEW(J,12)=OLD(J,14)
          NEW(J,13)=OLD(J,3)
          NEW(J,14)=OLD(J,15)
          NEW(J,15)=OLD(J,16)
          NEW(J,16)=OLD(J,4)
        ELSE
           STOP
        END IF
      ELSE IF (BASE_INFO%HEX_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==2) THEN
          !3D HEX LINEAR
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
          NEW(J,5)=OLD(J,5)
          NEW(J,6)=OLD(J,6)
          NEW(J,7)=OLD(J,7)
          NEW(J,8)=OLD(J,8)
        ELSE IF(BASE_INFO%B(I)%n==3) THEN
          !3D HEX QUADR
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,9)
          NEW(J,3)=OLD(J,2)
          NEW(J,4)=OLD(J,10)
          NEW(J,5)=OLD(J,11)
          NEW(J,6)=OLD(J,12)
          NEW(J,7)=OLD(J,3)
          NEW(J,8)=OLD(J,13)
          NEW(J,9)=OLD(J,4)
          NEW(J,10)=OLD(J,14)
          NEW(J,11)=OLD(J,15)
          NEW(J,12)=OLD(J,16)
          NEW(J,13)=OLD(J,17)
          NEW(J,14)=OLD(J,18)
          NEW(J,15)=OLD(J,19)
          NEW(J,16)=OLD(J,20)
          NEW(J,17)=OLD(J,21)
          NEW(J,18)=OLD(J,22)
          NEW(J,19)=OLD(J,5)
          NEW(J,20)=OLD(J,23)
          NEW(J,21)=OLD(J,6)
          NEW(J,22)=OLD(J,24)
          NEW(J,23)=OLD(J,25)
          NEW(J,24)=OLD(J,26)
          NEW(J,25)=OLD(J,7)
          NEW(J,26)=OLD(J,27)
          NEW(J,27)=OLD(J,8)
        ELSE IF(BASE_INFO%B(I)%n==4) THEN
          !3D HEX CUB
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,9)
          NEW(J,3)=OLD(J,10)
          NEW(J,4)=OLD(J,2)
          NEW(J,5)=OLD(J,11)
          NEW(J,6)=OLD(J,12)
          NEW(J,7)=OLD(J,13)
          NEW(J,8)=OLD(J,14)
          NEW(J,9)=OLD(J,15)
          NEW(J,10)=OLD(J,16)
          NEW(J,11)=OLD(J,17)
          NEW(J,12)=OLD(J,18)
          NEW(J,13)=OLD(J,3)
          NEW(J,14)=OLD(J,19)
          NEW(J,15)=OLD(J,20)
          NEW(J,16)=OLD(J,4)
          NEW(J,17)=OLD(J,21)
          NEW(J,18)=OLD(J,22)
          NEW(J,19)=OLD(J,23)
          NEW(J,20)=OLD(J,24)
          NEW(J,21)=OLD(J,25)
          NEW(J,22)=OLD(J,26)
          NEW(J,23)=OLD(J,27)
          NEW(J,24)=OLD(J,28)
          NEW(J,25)=OLD(J,29)
          NEW(J,26)=OLD(J,30)
          NEW(J,27)=OLD(J,31)
          NEW(J,28)=OLD(J,32)
          NEW(J,29)=OLD(J,33)
          NEW(J,30)=OLD(J,34)
          NEW(J,31)=OLD(J,35)
          NEW(J,32)=OLD(J,36)
          NEW(J,33)=OLD(J,37)
          NEW(J,34)=OLD(J,38)
          NEW(J,35)=OLD(J,39)
          NEW(J,36)=OLD(J,40)
          NEW(J,37)=OLD(J,41)
          NEW(J,38)=OLD(J,42)
          NEW(J,39)=OLD(J,43)
          NEW(J,40)=OLD(J,44)
          NEW(J,41)=OLD(J,45)
          NEW(J,42)=OLD(J,46)
          NEW(J,43)=OLD(J,47)
          NEW(J,44)=OLD(J,48)
          NEW(J,45)=OLD(J,49)
          NEW(J,46)=OLD(J,50)
          NEW(J,47)=OLD(J,51)
          NEW(J,48)=OLD(J,52)
          NEW(J,49)=OLD(J,5)
          NEW(J,50)=OLD(J,53)
          NEW(J,51)=OLD(J,54)
          NEW(J,52)=OLD(J,6)
          NEW(J,53)=OLD(J,55)
          NEW(J,54)=OLD(J,56)
          NEW(J,55)=OLD(J,57)
          NEW(J,56)=OLD(J,58)
          NEW(J,57)=OLD(J,59)
          NEW(J,58)=OLD(J,60)
          NEW(J,59)=OLD(J,61)
          NEW(J,60)=OLD(J,62)
          NEW(J,61)=OLD(J,7)
          NEW(J,62)=OLD(J,63)
          NEW(J,63)=OLD(J,64)
          NEW(J,64)=OLD(J,8)        
        ELSE
          STOP
        END IF
      ELSE IF (BASE_INFO%TRI_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==3) THEN
          !2D TET LINEAR
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
        ELSE IF(BASE_INFO%B(I)%n==6) THEN
          !2D TET QUAD
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
          NEW(J,5)=OLD(J,6)
          NEW(J,6)=OLD(J,5)
        ELSE IF(BASE_INFO%B(I)%n==10) THEN
          !2D TET CUB
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
          NEW(J,5)=OLD(J,5)
          NEW(J,6)=OLD(J,8)
          NEW(J,7)=OLD(J,9)
          NEW(J,8)=OLD(J,7)
          NEW(J,9)=OLD(J,6)
          NEW(J,10)=OLD(J,10)
        ELSE
          STOP
        END IF
      ELSE IF (BASE_INFO%TET_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==4) THEN
          !3D TET LINEAR
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==10) THEN
          !3D TET QUAD
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
          NEW(J,5)=OLD(J,5)
          NEW(J,6)=OLD(J,6)
          NEW(J,7)=OLD(J,7)
          NEW(J,8)=OLD(J,8)
          NEW(J,9)=OLD(J,10)
          NEW(J,10)=OLD(J,9)
        ELSE IF(BASE_INFO%B(I)%n==20) THEN
          !3D TET CUB
          NEW(J,1)=OLD(J,1)
          NEW(J,2)=OLD(J,2)
          NEW(J,3)=OLD(J,3)
          NEW(J,4)=OLD(J,4)
          NEW(J,5)=OLD(J,5)
          NEW(J,6)=OLD(J,6)
          NEW(J,7)=OLD(J,7)
          NEW(J,8)=OLD(J,8)
          NEW(J,9)=OLD(J,9)
          NEW(J,10)=OLD(J,10)
          NEW(J,11)=OLD(J,11)
          NEW(J,12)=OLD(J,12)
          NEW(J,13)=OLD(J,15)
          NEW(J,14)=OLD(J,16)
          NEW(J,15)=OLD(J,13)
          NEW(J,16)=OLD(J,14)
          NEW(J,17)=OLD(J,17)
          NEW(J,18)=OLD(J,18)
          NEW(J,19)=OLD(J,19)
          NEW(J,20)=OLD(J,20)
        ELSE
          STOP
        END IF
      ELSE
        STOP
      END IF
    END DO
!     RETURN
! 999 CALL ERRORS("FLUID_MECHANICS_IO_ORDER_NUMBERING",ERR,ERROR)    
!     RETURN
  END SUBROUTINE FLUID_MECHANICS_IO_ORDER_NUMBERING

  ! OK
  !================================================================================================================================
  !

  !> Combines nodes defined separately in cmheart.
  SUBROUTINE FLUID_MECHANICS_IO_MAKE_UNIQUE

!     INTEGER(INTG) :: ERR
!     TYPE(VARYING_STRING):: ERROR
  
    MESH_INFO(1)%T=MESH_INFO(1)%T
    MESH_INFO(2)%T=MESH_INFO(2)%T+ArrayOfNodesDefined(1)
    MESH_INFO(3)%T=MESH_INFO(3)%T+ArrayOfNodesDefined(1)+ArrayOfNodesDefined(2)
    
    IF(ArrayOfNodesDefined(1)==ArrayOfNodesDefined(2)) THEN
      ! copy all node numbers from 2 -> 1
      MESH_INFO(2)%T(:,:)= MESH_INFO(1)%T(:,:)
    ELSE
      IF (BASE_INFO%TRI_BASIS== 1) THEN
        MESH_INFO(2)%T(:,1:3)=MESH_INFO(1)%T(:,1:3)
      ELSE IF (BASE_INFO%TET_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%QUAD_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%HEX_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:8)=MESH_INFO(1)%T(:,1:8)
      ELSE
        STOP
      END IF
    END IF
    
    IF(ArrayOfNodesDefined(1)==ArrayOfNodesDefined(3)) THEN
      ! copy all node numbers from 2 -> 1
      MESH_INFO(3)%T(:,:)=MESH_INFO(1)%T(:,:)
    ELSE IF(ArrayOfNodesDefined(2)==ArrayOfNodesDefined(3)) THEN
      MESH_INFO(3)%T(:,:)=MESH_INFO(2)%T(:,:)
    ELSE
      IF (BASE_INFO%TRI_BASIS== 1) THEN
        MESH_INFO(3)%T(:,1:3)=MESH_INFO(1)%T(:,1:3)
      ELSE IF (BASE_INFO%TET_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%QUAD_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%HEX_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:8)=MESH_INFO(1)%T(:,1:8)
      ELSE
        STOP
      END IF
    END IF

!     RETURN
! 999 CALL ERRORS("FLUID_MECHANICS_IO_MAKE_UNIQUE",ERR,ERROR)    
!     RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_MAKE_UNIQUE
    
  !   OK
  !================================================================================================================================
  !

  !> Executes cmheart node reading process.
  SUBROUTINE FLUID_MECHANICS_IO_READ_NODES

    INTEGER(INTG):: I,J
    INTEGER(INTG):: a,b
!     INTEGER(INTG) :: ERR
!     TYPE(VARYING_STRING):: ERROR
    REAL(DP) :: TEMP(3)
!   REAL(DP),DIMENSION(832,3):: sebo_test_array
!   sebo_test_array=0.0

    READ(42,*) NAMz
    OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
    READ(1,*) ArrayOfNodesDefined(1:3)

    TotalNumberOfNodes=ArrayOfNodesDefined(1)+ArrayOfNodesDefined(2)+ArrayOfNodesDefined(3)
! ALLOCATE AND READ MESH NODE INFORMATION
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Reading Nodes...",ERR,ERROR,*999)
    WRITE(*,*)'Reading Nodes...'
    DO I=1,3
      MESH_INFO(I)%Lx=ArrayOfNodesDefined(I)
      ALLOCATE(MESH_INFO(I)%X(MESH_INFO(I)%Lx,3),STAT=ALLOC_ERROR)
      DO J = 1,MESH_INFO(I)%Lx
        READ(1,*,END=35) TEMP(1:3)
        MESH_INFO(I)%X(J,1:3)=TEMP(1:3)
!	WRITE(*,*) MESH_INFO(I)%X(J,1:3)
!        READ(1,*,END=35) sebo_test_array(J,1:3)
!	sebo_test_array(J,1:3)=(/1,2,3/)
      END DO
    END DO
    CLOSE(1)

    IF(ALLOCATED(OPENCMISS_NODE_COORD)) DEALLOCATE(OPENCMISS_NODE_COORD)
    IF(.NOT.ALLOCATED(OPENCMISS_NODE_COORD)) ALLOCATE(OPENCMISS_NODE_COORD(TotalNumberOfNodes,3),STAT=ALLOC_ERROR)
    a=1
    b=0
    DO I=1,3
      a=b+1
      b=b+ArrayOfNodesDefined(I)
      OPENCMISS_NODE_COORD(a:b,1:3)=MESH_INFO(I)%X(1:ArrayOfNodesDefined(I),1:3)
    END DO
    RETURN

35  PRINT *, 'FAILS'
    STOP
    IF(ALLOC_ERROR.NE.0) THEN
      STOP 'Error during allocation'
    END IF

    RETURN
! 999 CALL ERRORS("FLUID_MECHANICS_IO_READ_NODES",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_NODES

  ! OK
  !================================================================================================================================
  !
  !> Reads boundary conditions from a file
  SUBROUTINE FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION(SOLVER_TYPE,BOUNDARY_VALUES,BOUNDARY_NODES, &
    & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION,OPTION,ITERATION)
    INTEGER(INTG):: SOLVER_TYPE,I,NUMBER_OF_TIME_STEPS,OPTION
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    INTEGER(INTG), POINTER :: BOUNDARY_NODES(:)
    INTEGER(INTG):: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION,NUM_BC_NODES,ITERATION

    CHARACTER(34) :: INPUT_FILE

    INTEGER(INTG):: ENDI

    OPTION=1
    ENDI=1
    NUMBER_OF_TIME_STEPS=1
    NUMBER_OF_DIMENSIONS=42

    IF(SOLVER_TYPE==1) THEN !LINEAR
      IF(BOUNDARY_CONDITION==1)THEN !SCALAR BOUNDARY CONDITION - SET THE NUMBER OF DIMENSIONS TO BE NUMBER OF DEPENDENT FIELD COMPONENTS (this will usually be one)!  
       !ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS 
!       IF(J<10) THEN
          WRITE(INPUT_FILE,'("./input/BC/BC_VALUES_",I0,".dat")') ITERATION
          !WRITE (*,*) INPUT_FILE
!        ELSE IF(J<100) THEN
!          WRITE(INPUT_FILE,'("./input/BC/BC_VALUES_",F10,".dat")') J
!        ENDIF
        OPEN(UNIT=1, FILE=INPUT_FILE,STATUS='unknown')
          READ(1,*) NUM_BC_NODES
          ALLOCATE(BOUNDARY_NODES(NUM_BC_NODES))
          ALLOCATE(BOUNDARY_VALUES(NUM_BC_NODES))
        DO I=1,NUM_BC_NODES
          READ(1,*) BOUNDARY_NODES(I)
        ENDDO
        DO I=1,NUM_BC_NODES
          READ(1,*) BOUNDARY_VALUES(I)
        ENDDO
        BOUNDARY_VALUES=BOUNDARY_VALUES
        BOUNDARY_NODES=BOUNDARY_NODES
        !WRITE(*,*)'1! BOUNDARY_VALUES=BOUNDARY_VALUES'
        CLOSE(1)
      END IF
    ENDIF

    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION
  !
  !================================================================================================================================
  !
  !> Reads boundary conditions from a file
  SUBROUTINE FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_TYPE,BOUNDARY_VALUES,NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION, & 
    & OPTION,TIME_STEP,TIME,LENGTH_SCALE)

    INTEGER(INTG):: SOLVER_TYPE,I,OPTION, TIME_STEP
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    INTEGER(INTG):: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION,NUM_BC_NODES
    REAL(DP):: LENGTH_SCALE,TIME

    CHARACTER(34) :: INPUT_FILE

    !We assume a cardiac cycle to be 1s resolved with 20 time-steps which
    !gives us a time-step size of 0.05s

    INTEGER(INTG):: ENDI

    ENDI=42


    IF(SOLVER_TYPE==1) THEN !LINEAR
      IF(BOUNDARY_CONDITION==5)THEN !MOVED WALL
        IF(OPTION==0) THEN
          !do nothing (default)    
        ELSE IF(OPTION==1) THEN
          ENDI=SIZE(BOUNDARY_VALUES)
          IF(TIME_STEP<10) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
          ENDIF
          OPEN(UNIT=TIME_STEP,FILE=INPUT_FILE,STATUS='unknown') 
          DO I=1,ENDI
            READ(TIME_STEP,*) BOUNDARY_VALUES(I)
          ENDDO
          BOUNDARY_VALUES=BOUNDARY_VALUES*LENGTH_SCALE
          WRITE(*,*)'1! BOUNDARY_VALUES=BOUNDARY_VALUES'
          CLOSE(TIME_STEP)
        ELSEIF(OPTION==2) THEN
          ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS
          !U,X COMPONENT
          BOUNDARY_VALUES(1:ENDI)=0.0_DP                                       
          !V,Y COMPONENT
          BOUNDARY_VALUES(ENDI+1:ENDI+ENDI)=-0.05_DP*COS(2.0_DP*PI*1.0_DP/40.0_DP*TIME)
          !W,Z COMPONENT
          BOUNDARY_VALUES(ENDI+ENDI+1:ENDI+ENDI+ENDI)=0.0_DP
        ELSEIF(OPTION==3) THEN
          ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS
          !U,X COMPONENT
          BOUNDARY_VALUES(1:ENDI)=-1.0_DP
          !V,Y COMPONENT
          BOUNDARY_VALUES(ENDI+1:ENDI+ENDI)=-1.0_DP
          !W,Z COMPONENT
          BOUNDARY_VALUES(ENDI+ENDI+1:ENDI+ENDI+ENDI)=-1.0_DP
        ELSEIF(OPTION==4) THEN
          ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS
          !U,X COMPONENT
          BOUNDARY_VALUES(1:ENDI)=0.0_DP                                       
          !V,Y COMPONENT
          BOUNDARY_VALUES(ENDI+1:ENDI+ENDI)=0.0_DP
          !W,Z COMPONENT
          BOUNDARY_VALUES(ENDI+ENDI+1:ENDI+ENDI+ENDI)=-0.1_DP*COS(2.0_DP*PI*1.0_DP/40.0_DP*TIME)
        ELSEIF(OPTION==69) THEN
          ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS
          !U,X COMPONENT
          BOUNDARY_VALUES(1:ENDI)=0.0_DP                                       
          !V,Y COMPONENT
          BOUNDARY_VALUES(ENDI+1:ENDI+ENDI)=0.0_DP
          !W,Z COMPONENT
          BOUNDARY_VALUES(ENDI+ENDI+1:ENDI+ENDI+ENDI)=-100.0_DP*TIME
        ELSE
          STOP 'Error during boundary input'
        ENDIF
      ELSEIF(BOUNDARY_CONDITION==1)THEN !SCALAR BOUNDARY CONDITION - SET THE NUMBER OF DIMENSIONS TO BE NUMBER OF DEPENDENT FIELD COMPONENTS (this will usually be one)!  
       !ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS 
!       IF(J<10) THEN
          WRITE(INPUT_FILE,'("./input/BC/BC_VALUES_",F10.7,".dat")') TIME
          !WRITE (*,*) INPUT_FILE
!        ELSE IF(J<100) THEN
!          WRITE(INPUT_FILE,'("./input/BC/BC_VALUES_",F10,".dat")') J
!        ENDIF
        OPEN(UNIT=1, FILE=INPUT_FILE,STATUS='unknown')
          READ(1,*) NUM_BC_NODES
          ALLOCATE(BOUNDARY_VALUES(NUM_BC_NODES))
        DO I=1,NUM_BC_NODES
          READ(1,*) BOUNDARY_VALUES(I)
        ENDDO
        BOUNDARY_VALUES=BOUNDARY_VALUES
        WRITE(*,*)'1! BOUNDARY_VALUES=BOUNDARY_VALUES'
        CLOSE(1)
      END IF
    ELSE IF(SOLVER_TYPE==2) THEN !NONLINEAR
      IF(BOUNDARY_CONDITION==2)THEN !FIXED INLET
        IF(OPTION==69) THEN
          ENDI=SIZE(BOUNDARY_VALUES)/NUMBER_OF_DIMENSIONS
          !U,X COMPONENT
          BOUNDARY_VALUES(1:ENDI)=0.0_DP                                       
          !V,Y COMPONENT
          BOUNDARY_VALUES(ENDI+1:ENDI+ENDI)=0.0_DP
          !W,Z COMPONENT
          BOUNDARY_VALUES(ENDI+ENDI+1:ENDI+ENDI+ENDI)=1000.0_DP*TIME
        ELSE
          STOP 'Error during boundary input'
        ENDIF
      END IF
    ELSE IF(SOLVER_TYPE==3) THEN
      IF(BOUNDARY_CONDITION==2)THEN !FIXED INLET
     !do nothing
      END IF
    ENDIF
!         WRITE(*,*)'TEST_OUTPUT',BOUNDARY_VALUES(ENDI+1),PI,TIME
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS 

  ! OK
  !================================================================================================================================
  !

  !> Reads input data from a file
  SUBROUTINE FLUID_MECHANICS_IO_READ_DATA(SOLVER_TYPE,INPUT_VALUES,NUMBER_OF_DIMENSIONS,INPUT_TYPE, & 
    & INPUT_OPTION,TIME_STEP,LENGTH_SCALE)

    INTEGER(INTG):: SOLVER_TYPE,I,INPUT_OPTION,CHECK
    INTEGER(INTG) :: ERR
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    TYPE(VARYING_STRING):: ERROR
    REAL(DP), POINTER :: INPUT_VALUES(:)
    INTEGER(INTG):: INPUT_TYPE
    REAL(DP) :: LENGTH_SCALE

    INTEGER(INTG):: ENDI,TIME_STEP
    CHARACTER(35) :: INPUT_FILE
    CHARACTER(29) :: UVEL_FILE

    I=SOLVER_TYPE
    NUMBER_OF_DIMENSIONS=42  
    ENDI=42

  
!     IF(SOLVER_TYPE==1) THEN !LINEAR
      IF(INPUT_TYPE==1)THEN !POISSON VECTOR SOURCE TEMPORARY
        ENDI=SIZE(INPUT_VALUES)

! WRITE(*,*) "TIME_STEP", TIME_STEP

        IF(INPUT_OPTION==1) THEN
          IF(TIME_STEP<10) THEN
            WRITE(UVEL_FILE,'("./input/data/VEL_DATA_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(UVEL_FILE,'("./input/data/VEL_DATA_",I0,".dat")') TIME_STEP
          ENDIF
        ELSE IF(INPUT_OPTION==2) THEN
          IF(TIME_STEP<=10) THEN
            WRITE(UVEL_FILE,'("./input/data/VEL_DATA_0",I0,".dat")') TIME_STEP-1
          ELSE IF(TIME_STEP<100) THEN
            WRITE(UVEL_FILE,'("./input/data/VEL_DATA_",I0,".dat")') TIME_STEP-1
          ENDIF
        ELSE IF(INPUT_OPTION==3) THEN
          IF(TIME_STEP<10) THEN
            WRITE(UVEL_FILE,'("./input/data/ORI_DATA_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(UVEL_FILE,'("./input/data/ORI_DATA_",I0,".dat")') TIME_STEP
          ENDIF
        ELSE IF(INPUT_OPTION==4) THEN
          IF(TIME_STEP<10) THEN
            WRITE(UVEL_FILE,'("./input/data/U_DATA_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(UVEL_FILE,'("./input/data/U_DATA_",I0,".dat")') TIME_STEP
          ENDIF
        ELSE IF(INPUT_OPTION==5) THEN
          IF(TIME_STEP<10) THEN
            WRITE(UVEL_FILE,'("./input/data/V_DATA_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(UVEL_FILE,'("./input/data/V_DATA_",I0,".dat")') TIME_STEP
          ENDIF
        ELSE IF(INPUT_OPTION==6) THEN
          IF(TIME_STEP<10) THEN
            WRITE(UVEL_FILE,'("./input/data/W_DATA_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(UVEL_FILE,'("./input/data/W_DATA_",I0,".dat")') TIME_STEP
          ENDIF
        ENDIF
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,UVEL_FILE,ERR,ERROR,*999)
          OPEN(UNIT=42, FILE=UVEL_FILE,STATUS='unknown') 
          READ(42,*) CHECK
          IF(CHECK/=ENDI) THEN
            STOP 'Error during data input - probably wrong Lagrangian/Hermite input file!'
          ENDIF
          DO I=1,ENDI
            READ(42,*) INPUT_VALUES(I)
          ENDDO
          CLOSE(42)

      ELSE IF(INPUT_TYPE==42) THEN
! do nothing for now
        IF(INPUT_OPTION==0) THEN
          !do nothing (default)    
        ELSE IF(INPUT_OPTION==1) THEN
          ENDI=SIZE(INPUT_VALUES)
          IF(TIME_STEP<=10) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_00",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<1000) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
          ENDIF
          OPEN(UNIT=TIME_STEP, FILE=INPUT_FILE,STATUS='unknown') 
          DO I=1,ENDI
            READ(TIME_STEP,*) INPUT_VALUES(I)
          ENDDO
! ! TESTETSTEST
          INPUT_VALUES=INPUT_VALUES/LENGTH_SCALE
          WRITE(*,*)'1! INPUT_VALUES=INPUT_VALUES/LENGTH_SCALE'
          CLOSE(TIME_STEP)
        ELSE IF(INPUT_OPTION==2) THEN ! For Darcy, invoke the length scale (consistent with reading in the geometry data)
          ENDI=SIZE(INPUT_VALUES)
          IF(TIME_STEP<=10) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_00",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<1000) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
          ENDIF
          OPEN(UNIT=TIME_STEP, FILE=INPUT_FILE,STATUS='unknown') 
          DO I=1,ENDI
            READ(TIME_STEP,*) INPUT_VALUES(I)
            INPUT_VALUES(I) = LENGTH_SCALE * INPUT_VALUES(I)
          ENDDO
          CLOSE(TIME_STEP)
        ELSE IF(INPUT_OPTION==3) THEN
          ENDI=SIZE(INPUT_VALUES)
          IF(TIME_STEP<=10) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_00",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<100) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<500) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
          ELSE IF(TIME_STEP<1000) THEN
            WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') 1000-TIME_STEP
          ENDIF
          OPEN(UNIT=TIME_STEP, FILE=INPUT_FILE,STATUS='unknown') 
          DO I=1,ENDI
            READ(TIME_STEP,*) INPUT_VALUES(I)
            INPUT_VALUES(I) = LENGTH_SCALE * INPUT_VALUES(I)
          ENDDO
          CLOSE(TIME_STEP)
        ELSE
          STOP 'Error during data input'
        ENDIF
      ENDIF
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_IO_READ_DATA",ERR,ERROR)    
    RETURN
  END SUBROUTINE FLUID_MECHANICS_IO_READ_DATA

  ! OK
  !================================================================================================================================
  !

  !> Execute cmheart element reading process.
  SUBROUTINE FLUID_MECHANICS_IO_READ_ELEMENTS
   
    INTEGER(INTG):: I,J
!     INTEGER(INTG) :: ERR
!     TYPE(VARYING_STRING):: ERROR
    INTEGER(INTG) :: TEMP(200)

    READ(42,*) NAMz
    CLOSE(42)
    OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
    READ(1,*) NumberOfElementsDefined(1:3)

! ALLOCATE AND READ MESH ELEMENT INFORMATION
    TotalNumberOfNodes=ArrayOfNodesDefined(1)+ArrayOfNodesDefined(2)+ArrayOfNodesDefined(3)
! ! !     CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Reading Elements...",ERR,ERROR,*999)
    WRITE(*,*)'Reading Elements...'            
    DO I=1,3
      MESH_INFO(I)%Lt=NumberOfElementsDefined(I)
      ALLOCATE(MESH_INFO(I)%T(MESH_INFO(I)%Lt,NumberOfNodesPerElement(I)),STAT=ALLOC_ERROR)
        DO J = 1,MESH_INFO(I)%Lt
          READ(1,*,END=30) TEMP(1:NumberOfNodesPerElement(I))
          MESH_INFO(I)%T(J,1:NumberOfNodesPerElement(I))=TEMP(1:NumberOfNodesPerElement(I))
        END DO
    END DO
    CLOSE(1)
    RETURN

 30 PRINT *, 'FAILS'
    STOP

    IF(ALLOC_ERROR.NE.0) THEN
      STOP 'Error during allocation'
    END IF

    RETURN
! 999 CALL ERRORS("FLUID_MECHANICS_IO_READ_ELEMENTS",ERR,ERROR)    
    RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_READ_ELEMENTS

  ! OK
  !================================================================================================================================
  !

  !> Writes nice information on screen.
  SUBROUTINE FLUID_MECHANICS_IO_PRINT_ON_SCREEN

    INTEGER(INTG):: I
!     INTEGER(INTG) :: ERR
!     TYPE(VARYING_STRING):: ERROR

    DO I = 1,TotalNumberOfNodes
      WRITE(*,'("Node ",(I0,4x),1000( F5.3,2x ))')I,OPENCMISS_NODE_COORD(I,1:3)
    END DO

  ! where are the element nodes stored -> 3 MATRICES

    WRITE(*,*)
    WRITE(*,*)
    DO I = 1,NumberOfElementsDefined(1)
      WRITE(*,'("M-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_M(I,1:NumberOfNodesPerElement(1))
    END DO
    WRITE(*,*)
    DO I = 1,NumberOfElementsDefined(2)
      WRITE(*,'("V-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_V(I,1:NumberOfNodesPerElement(2))
    END DO
    WRITE(*,*)
    DO I = 1,NumberOfElementsDefined(3)
      WRITE(*,'("P-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_P(I,1:NumberOfNodesPerElement(3))
    END DO
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)

!     RETURN
! 999 CALL ERRORS("FLUID_MECHANICS_IO_PRINT_ON_SCREEN",ERR,ERROR)    
!     RETURN

  END SUBROUTINE FLUID_MECHANICS_IO_PRINT_ON_SCREEN

  ! OK
  !================================================================================================================================
  !

  !================================================================================================================================
  ! 
  ! Below are the Darcy routines
  ! chrm, 20.08.09
  !
  !================================================================================================================================

  ! OK
  !================================================================================================================================
  !

  !> Read parameters for Darcy problem
  SUBROUTINE FLUID_MECHANICS_IO_READ_DARCY_PARAMS

    IMPLICIT NONE

    CHARACTER*90 DARCY_PARAM_FILE

    DARCY_PARAM_FILE='./input/Darcy_parameters.inp'

    WRITE(*,*)'Reading Darcy parameters.'

    OPEN(UNIT=37,FILE=DARCY_PARAM_FILE,STATUS='old',action='read') ! Read base file for initial parameters


    DO WHILE (0 < 1)
      READ(37,*,END=50) IN_CHAR
      IF (INDEX(IN_CHAR,'TESTCASE:') == 1)      READ(37,*) DARCY%TESTCASE

      IF (INDEX(IN_CHAR,'STAB:') == 1)          READ(37,*) DARCY%STAB
      IF (INDEX(IN_CHAR,'ANALYTIC:') == 1)      READ(37,*) DARCY%ANALYTIC
      IF (INDEX(IN_CHAR,'DEBUG:') == 1)         READ(37,*) DARCY%DEBUG

      IF (INDEX(IN_CHAR,'LENGTH:') == 1)        READ(37,*) DARCY%LENGTH

      IF (INDEX(IN_CHAR,'GEOM_TOL:') == 1)      READ(37,*) DARCY%GEOM_TOL
      IF (INDEX(IN_CHAR,'X1:') == 1)            READ(37,*) DARCY%X1
      IF (INDEX(IN_CHAR,'X2:') == 1)            READ(37,*) DARCY%X2
      IF (INDEX(IN_CHAR,'Y1:') == 1)            READ(37,*) DARCY%Y1
      IF (INDEX(IN_CHAR,'Y2:') == 1)            READ(37,*) DARCY%Y2
      IF (INDEX(IN_CHAR,'Z1:') == 1)            READ(37,*) DARCY%Z1
      IF (INDEX(IN_CHAR,'Z2:') == 1)            READ(37,*) DARCY%Z2

      IF (INDEX(IN_CHAR,'PERM:') == 1)          READ(37,*) DARCY%PERM
      IF (INDEX(IN_CHAR,'VIS:') == 1)           READ(37,*) DARCY%VIS
      IF (INDEX(IN_CHAR,'P_SINK:') == 1)        READ(37,*) DARCY%P_SINK

      IF (INDEX(IN_CHAR,'BC_NUMBER_OF_WALL_NODES:') == 1) READ(37,*) DARCY%BC_NUMBER_OF_WALL_NODES
      IF (INDEX(IN_CHAR,'NUMBER_OF_BCS:') == 1) READ(37,*) DARCY%NUMBER_OF_BCS
    END DO

50 CLOSE(37)

    IF( DARCY%DEBUG ) THEN
      write(*,*)'Read Darcy parameters from the following file = ',DARCY_PARAM_FILE
      write(*,*)'Press ENTER to continue.'
      read(*,*)
    END IF

    IF( ABS(DARCY%VIS) > 1.0E-14 ) THEN
      DARCY%PERM_OVER_VIS = DARCY%PERM / DARCY%VIS
    ELSE
      WRITE(*,*)'Darcy_parameters: VIS cannot be machine zero.'
      STOP
    END IF

    DARCY%max_node_spacing = 0.05_DP

  END SUBROUTINE FLUID_MECHANICS_IO_READ_DARCY_PARAMS

  ! OK
  !================================================================================================================================
  !

  !> Get analytic solution for Darcy problem
  SUBROUTINE FLUID_MECHANICS_IO_DARCY_GET_ANALYTIC

    INTEGER(INTG):: I
    REAL(DP):: COORD_X, COORD_Y, COORD_Z, ARG_X, ARG_Y, ARG_Z
    REAL(DP):: FACT

    IF( DARCY%TESTCASE == 1 ) THEN
      FACT = DARCY%PERM_OVER_VIS
    ELSE IF( DARCY%TESTCASE == 2 ) THEN
      FACT = 2.0_DP * PI * DARCY%PERM_OVER_VIS / DARCY%LENGTH
!       FACT = 1.0_DP
    ELSE IF( DARCY%TESTCASE == 3 ) THEN
      FACT = - DARCY%PERM_OVER_VIS * 2.0_DP * PI / DARCY%LENGTH
    END IF

    IF( NumberOfDimensions==2 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        COORD_X = NodeXValue(I)
        COORD_Y = NodeYValue(I)

        IF( DARCY%TESTCASE == 1 ) THEN
          NodeUValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X + 2.0_DP * COORD_Y )
          NodeVValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X - 2.0_DP * COORD_Y )
          NodePValue_analytic(I) = COORD_X * COORD_X + 2.0_DP * COORD_X * COORD_Y - COORD_Y * COORD_Y
        ELSE IF( DARCY%TESTCASE == 2 ) THEN
          ARG_X = 2.0_DP * PI * COORD_X / DARCY%LENGTH
          ARG_Y = 2.0_DP * PI * COORD_Y / DARCY%LENGTH
          NodeUValue_analytic(I) = - FACT * COS( ARG_X ) * SIN( ARG_Y ) 
          NodeVValue_analytic(I) = - FACT * SIN( ARG_X ) * COS( ARG_Y ) 
          NodePValue_analytic(I) =          SIN( ARG_X ) * SIN( ARG_Y )
        ELSE IF( DARCY%TESTCASE == 3 ) THEN
          ARG_X = 2.0_DP * PI * COORD_X / DARCY%LENGTH
          ARG_Y = 2.0_DP * PI * COORD_Y / DARCY%LENGTH
          NodeUValue_analytic(I) = FACT  * ( 9.0_DP * COS( ARG_X ) )
          NodeVValue_analytic(I) = FACT  * ( 1.0_DP * SIN( ARG_Y ) )
          NodePValue_analytic(I) =           9.0_DP * SIN( ARG_X ) - 1.0_DP * COS( ARG_Y ) + DARCY%P_SINK
        END IF
      END DO
    ELSE IF( NumberOfDimensions==3 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        COORD_X = NodeXValue(I)
        COORD_Y = NodeYValue(I)
        COORD_Z = NodeZValue(I)

        IF( DARCY%TESTCASE == 1 ) THEN
          NodeUValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X + 2.0_DP * COORD_Y + COORD_Z )
          NodeVValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X - 2.0_DP * COORD_Y + COORD_Z )
          NodeWValue_analytic(I) = - FACT * ( 3.0_DP + COORD_X + COORD_Y )
          NodePValue_analytic(I) = COORD_X * COORD_X + 2.0_DP * COORD_X * COORD_Y - COORD_Y * COORD_Y + &
            & 3.0_DP * COORD_Z + COORD_Z * COORD_X  + COORD_Z * COORD_Y 
        ELSE IF( DARCY%TESTCASE == 2 ) THEN
          ARG_X = 2.0_DP * PI * COORD_X / DARCY%LENGTH
          ARG_Y = 2.0_DP * PI * COORD_Y / DARCY%LENGTH
          ARG_Z = 2.0_DP * PI * COORD_Z / DARCY%LENGTH
          NodeUValue_analytic(I) = - FACT * COS( ARG_X ) * SIN( ARG_Y )  * SIN( ARG_Z ) 
          NodeVValue_analytic(I) = - FACT * SIN( ARG_X ) * COS( ARG_Y )  * SIN( ARG_Z )  
          NodeWValue_analytic(I) = - FACT * SIN( ARG_X ) * SIN( ARG_Y )  * COS( ARG_Z )  
          NodePValue_analytic(I) =          SIN( ARG_X ) * SIN( ARG_Y )  * SIN( ARG_Z )  
        ELSE IF( DARCY%TESTCASE == 3 ) THEN
          ARG_X = 2.0_DP * PI * COORD_X / DARCY%LENGTH
          ARG_Y = 2.0_DP * PI * COORD_Y / DARCY%LENGTH
          ARG_Z = 2.0_DP * PI * COORD_Z / DARCY%LENGTH
          NodeUValue_analytic(I) = FACT  * ( 9.0_DP * COS( ARG_X ) )
          NodeVValue_analytic(I) = FACT  * ( 1.0_DP * SIN( ARG_Y ) )
          NodeWValue_analytic(I) = FACT  * (-3.0_DP * SIN( ARG_Z ) )
          NodePValue_analytic(I) =           9.0_DP * SIN( ARG_X ) - 1.0_DP * COS( ARG_Y ) &
            &                              + 3.0_DP * COS( ARG_Z ) + DARCY%P_SINK
        END IF
      END DO
    END IF

  END SUBROUTINE FLUID_MECHANICS_IO_DARCY_GET_ANALYTIC

  ! OK
  !================================================================================================================================
  !

  !> Evaluate local error for Darcy problem
  SUBROUTINE FLUID_MECHANICS_IO_DARCY_EVAL_ERROR

    INTEGER(INTG):: I

    IF( NumberOfDimensions==2 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        NodeUValue_error(I) = NodeUValue(I) - NodeUValue_analytic(I)
        NodeVValue_error(I) = NodeVValue(I) - NodeVValue_analytic(I)
        NodePValue_error(I) = NodePValue(I) - NodePValue_analytic(I)
      END DO
    ELSE IF( NumberOfDimensions==3 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        NodeUValue_error(I) = NodeUValue(I) - NodeUValue_analytic(I)
        NodeVValue_error(I) = NodeVValue(I) - NodeVValue_analytic(I)
        NodeWValue_error(I) = NodeWValue(I) - NodeWValue_analytic(I)
        NodePValue_error(I) = NodePValue(I) - NodePValue_analytic(I)
      END DO
    END IF

  END SUBROUTINE FLUID_MECHANICS_IO_DARCY_EVAL_ERROR

  ! OK
  !================================================================================================================================
  !

  !> Evaluate max error for Darcy problem
  SUBROUTINE FLUID_MECHANICS_IO_DARCY_EVAL_MAX_ERROR

    INTEGER(INTG):: I
    REAL(DP):: MaxNodeUValue_error, MaxNodeVValue_error, MaxNodeWValue_error, MaxNodePValue_error

    OPEN(UNIT=23, FILE='./output/conv.node',STATUS='unknown')

    MaxNodeUValue_error = 0.0
    MaxNodeVValue_error = 0.0
    MaxNodeWValue_error = 0.0
    MaxNodePValue_error = 0.0

    IF( NumberOfDimensions==2 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        IF( abs(mod( ((NodeXValue(I)-DARCY%X1) / DARCY%max_node_spacing), 1.0_DP)) < DARCY%GEOM_TOL ) THEN
          IF( abs(mod( ((NodeYValue(I)-DARCY%Y1) / DARCY%max_node_spacing), 1.0_DP)) < DARCY%GEOM_TOL ) THEN

              WRITE(23,'("    ", es25.16 )')NodeXValue(I)
              WRITE(23,'("    ", es25.16 )')NodeYValue(I)

              WRITE(23,'("    ", es25.16 )')NodeUValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodeVValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodePValue_error(I)

              WRITE(23,*) ' '

              IF( abs(NodeUValue_error(I)) > MaxNodeUValue_error ) MaxNodeUValue_error = abs(NodeUValue_error(I))
              IF( abs(NodeVValue_error(I)) > MaxNodeVValue_error ) MaxNodeVValue_error = abs(NodeVValue_error(I))
              IF( abs(NodePValue_error(I)) > MaxNodePValue_error ) MaxNodePValue_error = abs(NodePValue_error(I))

          END IF
        END IF
      END DO
      WRITE(23,'("    MaxNodeUValue_error = ", es25.16 )')MaxNodeUValue_error
      WRITE(23,'("    MaxNodeVValue_error = ", es25.16 )')MaxNodeVValue_error
      WRITE(23,'("    MaxNodePValue_error = ", es25.16 )')MaxNodePValue_error
      WRITE(23,*) ' '
    ELSE IF( NumberOfDimensions==3 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        IF( abs(mod( ((NodeXValue(I)-DARCY%X1) / DARCY%max_node_spacing), 1.0_DP)) < DARCY%GEOM_TOL ) THEN
          IF( abs(mod( ((NodeYValue(I)-DARCY%Y1) / DARCY%max_node_spacing), 1.0_DP)) < DARCY%GEOM_TOL ) THEN
            IF( abs(mod( ((NodeZValue(I)-DARCY%Z1) / DARCY%max_node_spacing), 1.0_DP)) < DARCY%GEOM_TOL ) THEN

              WRITE(23,'("    ", es25.16 )')NodeXValue(I)
              WRITE(23,'("    ", es25.16 )')NodeYValue(I)
              WRITE(23,'("    ", es25.16 )')NodeZValue(I)

              WRITE(23,'("    ", es25.16 )')NodeUValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodeVValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodeWValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodePValue_error(I)

              WRITE(23,*) ' '

              IF( abs(NodeUValue_error(I)) > MaxNodeUValue_error ) MaxNodeUValue_error = abs(NodeUValue_error(I))
              IF( abs(NodeVValue_error(I)) > MaxNodeVValue_error ) MaxNodeVValue_error = abs(NodeVValue_error(I))
              IF( abs(NodeWValue_error(I)) > MaxNodeWValue_error ) MaxNodeWValue_error = abs(NodeWValue_error(I))
              IF( abs(NodePValue_error(I)) > MaxNodePValue_error ) MaxNodePValue_error = abs(NodePValue_error(I))

            END IF
          END IF
        END IF
      END DO
      WRITE(23,'("    MaxNodeUValue_error = ", es25.16 )')MaxNodeUValue_error
      WRITE(23,'("    MaxNodeVValue_error = ", es25.16 )')MaxNodeVValue_error
      IF( NumberOfDimensions==3 ) THEN
        WRITE(23,'("    MaxNodeWValue_error = ", es25.16 )')MaxNodeWValue_error
      END IF
      WRITE(23,'("    MaxNodePValue_error = ", es25.16 )')MaxNodePValue_error
      WRITE(23,*) ' '
    END IF

    CLOSE(23)

  END SUBROUTINE FLUID_MECHANICS_IO_DARCY_EVAL_MAX_ERROR

  ! OK
  !================================================================================================================================
  ! 

END MODULE FLUID_MECHANICS_IO_ROUTINES
