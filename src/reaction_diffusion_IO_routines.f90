!> \file
!> \author Vijay Rajagopal
!> \brief This module handles some mesh/parameter input routines and cmgui output routines for reaction diffusion
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

MODULE REACTION_DIFFUSION_IO_ROUTINES

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

  TYPE BOUNDARY_PARAMETERS
    INTEGER:: NUMBER_OF_BC
    INTEGER:: NUMBER_OF_BC1
    INTEGER:: NUMBER_OF_BC2
    INTEGER:: NUMBER_OF_BC3
    INTEGER:: NUMBER_OF_BC4
    INTEGER:: NUMBER_OF_BC5
    INTEGER, POINTER:: BC_NODE_NUMBER1(:)
    INTEGER, POINTER:: BC_NODE_NUMBER2(:)
    INTEGER, POINTER:: BC_NODE_NUMBER3(:)
    INTEGER, POINTER:: BC_NODE_NUMBER4(:)
    INTEGER, POINTER:: BC_NODE_NUMBER5(:)
  END TYPE BOUNDARY_PARAMETERS


  !Module variables

  TYPE (ARRAY_PROBLEM_BASE) BASE_INFO
  TYPE (ARRAY_MESH) MESH_INFO(3)
  TYPE (EXPORT_CONTAINER) TMP, TMP1
  TYPE (BOUNDARY_PARAMETERS) BC_TMP
  TYPE (COUPLING_PARAMETERS) CONNECT_TMP
  TYPE(FIELD_TYPE), POINTER :: FIELD, MATERIAL_FIELD, SOURCE_FIELD
  TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:),MATERIAL_INTERPOLATION_PARAMETERS(:)
  TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:),MATERIAL_INTERPOLATED_POINT(:)

  TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: SOURCE_INTERPOLATION_PARAMETERS(:)
  TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: SOURCE_INTERPOLATED_POINT(:)


  INTEGER(INTG), DIMENSION(:), ALLOCATABLE:: NodesPerElement
  INTEGER(INTG), DIMENSION(:,:), ALLOCATABLE::ElementNodes
  INTEGER(INTG):: NumberOfFields
  INTEGER(INTG):: NumberOfDimensions
  INTEGER(INTG):: ValueIndex
  INTEGER(INTG):: NumberOfVariableComponents
  INTEGER(INTG):: NumberOfMeshComponents
  INTEGER(INTG):: NumberOfMaterialComponents
  INTEGER(INTG):: NumberOfSourceComponents
  INTEGER(INTG):: NumberOfNodesDefined
  INTEGER(INTG):: NumberOfFieldComponent(4) !(3)
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

  LOGICAL :: OUTPUT_SOURCE

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
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeDIFFXValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeDIFFYValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeDIFFZValue
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeSTORAGEValue
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

  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeSourceValue1 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeSourceValue2 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeSourceValue3 
  REAL(DP), DIMENSION(:), ALLOCATABLE:: NodeSourceValue4 

  REAL(DP):: ScaleFactorsPerElementNodes(10,10)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE::OPENCMISS_NODE_COORD

  CHARACTER*2 NMs(99),KNOT
  CHARACTER*60 IN_CHAR
  CHARACTER*90 NIMZ
!   CHARACTER*30 NAMz
  CHARACTER*90 NAMz


  PUBLIC REACTION_DIFFUSION_IO_WRITE_CMGUI

CONTAINS

  ! OK
  !================================================================================================================================
  !

  !> Writes solution into cmgui formats exelem and exnode.
  SUBROUTINE REACTION_DIFFUSION_IO_WRITE_CMGUI(REGION, EQUATIONS_SET_GLOBAL_NUMBER, NAME, ERR, ERROR,*)

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

    CALL ENTERS("REACTION_DIFFUSION_IO_WRITE_CMGUI",ERR,ERROR,*999)

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
    IF (ALLOCATED(NodeDIFFXValue)) DEALLOCATE(NodeDIFFXValue) 
    IF (ALLOCATED(NodeDIFFYValue)) DEALLOCATE(NodeDIFFYValue) 
    IF (ALLOCATED(NodeDIFFZValue)) DEALLOCATE(NodeDIFFZValue) 
    IF (ALLOCATED(NodeSTORAGEValue)) DEALLOCATE(NodeSTORAGEValue) 
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

    IF (ALLOCATED(NodeSourceValue1)) DEALLOCATE(NodeSourceValue1) 
    IF (ALLOCATED(NodeSourceValue2)) DEALLOCATE(NodeSourceValue2) 
    IF (ALLOCATED(NodeSourceValue3)) DEALLOCATE(NodeSourceValue3) 
    IF (ALLOCATED(NodeSourceValue4)) DEALLOCATE(NodeSourceValue4) 

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
    WRITE(*,*) 'I AM IN REAC-DIFF-IO'
    EQUATIONS_SET => REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr

!---tob
!     FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE%VARIABLE_TYPE
!     ! '1' associated with linear matrix

    var_idx = 1
    FIELD_VAR_TYPE = FIELD_U_VARIABLE_TYPE
    parameter_set_idx = 1

!---toe

!    NumberOfFields=REGION%fields%number_of_fields
! Hack for ALE... to be removed later
    NumberOfFields=3
    NumberOfDimensions=REGION%coordinate_system%number_of_dimensions
    NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
      & variables(var_idx)%number_of_components

    NumberOfMaterialComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
      & variables(1)%number_of_components
    NumberOfElements=REGION%meshes%meshes(1)%ptr%number_of_elements
    NumberOfMeshComponents=REGION%meshes%meshes(1)%ptr%number_of_components

    IF(.NOT.ALLOCATED(NodesPerElement)) ALLOCATE(NodesPerElement(MAX(NumberOfMeshComponents,NumberOfElements)))
    
    IF(.NOT.ALLOCATED(NodesPerMeshComponent)) ALLOCATE(NodesPerMeshComponent(NumberOfMeshComponents))
    MaxNodesPerElement=0

    DO I=1,NumberOfMeshComponents
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(1)%basis%number_of_element_parameters
      NodesPerMeshComponent(I)=REGION%meshes%meshes(1)%ptr%topology(I)%ptr%nodes%number_of_nodes
    END DO

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
    IF(.NOT.ALLOCATED(NodeDIFFXValue)) ALLOCATE(NodeDIFFXValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeDIFFYValue)) ALLOCATE(NodeDIFFYValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeDIFFZValue)) ALLOCATE(NodeDIFFZValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSTORAGEValue)) ALLOCATE(NodeSTORAGEValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeMUValue)) ALLOCATE(NodeMUValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeLabelValue)) ALLOCATE(NodeLabelValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeRHOValue)) ALLOCATE(NodeRHOValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeA0Value)) ALLOCATE(NodeA0Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeH0Value)) ALLOCATE(NodeH0Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeEValue)) ALLOCATE(NodeEValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSIGMAValue)) ALLOCATE(NodeSIGMAValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeKappaValue)) ALLOCATE(NodeKappaValue(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,MaxNodesPerElement))
    IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,MaxNodesPerElement))

    IF(.NOT.ALLOCATED(NodePerm2Value)) ALLOCATE(NodePerm2Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm3Value)) ALLOCATE(NodePerm3Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm4Value)) ALLOCATE(NodePerm4Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm5Value)) ALLOCATE(NodePerm5Value(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodePerm6Value)) ALLOCATE(NodePerm6Value(NodesPerMeshComponent(1)))

    IF(.NOT.ALLOCATED(NodeSourceValue1)) ALLOCATE(NodeSourceValue1(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSourceValue2)) ALLOCATE(NodeSourceValue2(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSourceValue3)) ALLOCATE(NodeSourceValue3(NodesPerMeshComponent(1)))
    IF(.NOT.ALLOCATED(NodeSourceValue4)) ALLOCATE(NodeSourceValue4(NodesPerMeshComponent(1)))

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

    NULLIFY(MATERIAL_FIELD)
    NULLIFY(SOURCE_FIELD)

    FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field

    !material
    MATERIAL_FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field

    !source field
    OUTPUT_SOURCE = .FALSE.
    IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) & 
      & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
        & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
          SOURCE_FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field
          IF( ASSOCIATED(SOURCE_FIELD) ) OUTPUT_SOURCE = .TRUE.
    END IF

    NULLIFY(INTERPOLATION_PARAMETERS,MATERIAL_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATED_POINT,MATERIAL_INTERPOLATED_POINT,SOURCE_INTERPOLATED_POINT)

    CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)

    !material
    CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(MATERIAL_FIELD,MATERIAL_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    CALL FIELD_INTERPOLATED_POINTS_INITIALISE(MATERIAL_INTERPOLATION_PARAMETERS,MATERIAL_INTERPOLATED_POINT,ERR,ERROR,*999)

    !source field
    IF( OUTPUT_SOURCE ) THEN
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(SOURCE_FIELD,SOURCE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_INITIALISE(SOURCE_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATED_POINT,ERR,ERROR,*999)
      NumberOfSourceComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field% &
        & variables(1)%number_of_components

      NumberOfFields = NumberOfFields + 1
    END IF

    DO I=1,NumberOfElements

      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1) &
        & %ptr%topology%elements%elements(I)%basis%number_of_element_parameters

      DO J=1,NodesPerElement(I)
        ELEMENT_NUMBER=I

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


        !material
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & MATERIAL_INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,MATERIAL_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr, &
          & ERR,ERROR,*999)

        !source field
        IF( OUTPUT_SOURCE ) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & SOURCE_INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,SOURCE_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr, &
            & ERR,ERROR,*999)
        END IF


        NodeXValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)

        IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3)THEN
          NodeYValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))
        END IF
        IF(NumberOfDimensions==3)THEN
          NodeZValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
        END IF

        NodeUValue(K)=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
          & variables(var_idx)%parameter_sets%parameter_sets(parameter_set_idx)%ptr%parameters%cmiss%data_dp(K)

        MATERIAL_INTERPOLATION_TYPE=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
          & materials%materials_field%variables(1)%COMPONENTS(1)%INTERPOLATION_TYPE 

        IF(MATERIAL_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION)THEN

            NodeDIFFXValue(K)   =MATERIAL_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
            NodeDIFFYValue(K)  =MATERIAL_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1)
            NodeDIFFZValue(K)  =MATERIAL_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,1)
            NodeSTORAGEValue(K)  =MATERIAL_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(4,1)
            

          IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) & 
            & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
              & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
            !source field
            IF( OUTPUT_SOURCE ) THEN
              NodeSourceValue1(K)=SOURCE_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
            END IF
          END IF

        ELSE !default to FIELD_CONSTANT_INTERPOLATION
          NodeDIFFXValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%materials%materials_field% &
            & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1)
          IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3)THEN
            NodeDIFFYValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
              & materials%materials_field%variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(2)
          ENDIF
          IF(NumberOfDimensions==3)THEN
            NodeDIFFZValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
              & materials%materials_field%variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(3)
          ENDIF
          IF(NumberOfDimensions==1)THEN
            NodeSTORAGEValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
              & materials%materials_field%variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(2)
          ELSEIF(NumberOfDimensions==2)THEN
            NodeSTORAGEValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
              & materials%materials_field%variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(3)
          ELSE
            NodeSTORAGEValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr% &
              & materials%materials_field%variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(4)
          ENDIF




        END IF

      END DO 
    END DO

    DN=.TRUE.
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

    IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) & 
      & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
        & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
        IF( OUTPUT_SOURCE ) THEN
          NumberOfFieldComponent(4)=NumberOfSourceComponents
        END IF
    END IF

    DO I=1,NumberOfElements
      DO J=1,NodesPerElement(I)
        ElementNodes(I,J)=REGION%meshes%meshes(1)%ptr%topology(1)% &
          & ptr%elements%elements(I)%global_element_nodes(J)
        ElementNodesScales(I,J)=1.0000000000000000E+00
      END DO
    END DO

    CALL REACTION_DIFFUSION_IO_WRITE_NODES_CMGUI(NAME,EQUATIONS_SET)
    CALL REACTION_DIFFUSION_IO_WRITE_ELEMENTS_CMGUI(NAME)

    CALL EXITS("REACTION_DIFFUSION_IO_WRITE_CMGUI")
    RETURN
999 CALL ERRORS("REACTION_DIFFUSION_IO_WRITE_CMGUI",ERR,ERROR)    
    CALL EXITS("REACTION_DIFFUSION_IO_WRITE_CMGUI")
    RETURN

  END SUBROUTINE REACTION_DIFFUSION_IO_WRITE_CMGUI



  ! OK
  !================================================================================================================================
  !

  !> Executes nodes writing process.
  SUBROUTINE REACTION_DIFFUSION_IO_WRITE_NODES_CMGUI(NAME,EQUATIONS_SET)

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

    FILENAME="./output/"//NAME//".exnode"
    OPEN(UNIT=14, FILE=CHAR(FILENAME),STATUS='unknown')

! WRITING HEADER INFORMATION

    WRITE(14,*) 'Group name: 09ryrNtatsCell'

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

    WRITE(14,*) ' 2) dependent,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))

    DO I=1,NumberOfVariableComponents
      WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
      ValueIndex=ValueIndex+1
    END DO

    WRITE(14,*) ' 3) material,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfMaterialComponents))

    DO I=1,NumberOfMaterialComponents
      WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
      ValueIndex=ValueIndex+1
    END DO


    IF( OUTPUT_SOURCE ) THEN !Watch out that no numbering conflict occurs with Analytic: 4.)
      WRITE(14,*) ' 4) source,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfSourceComponents))

      DO I=1,NumberOfSourceComponents
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

      !write out materials
      
      WRITE(14,'("    ", es25.16 )')NodeDIFFXValue(I)
      IF(NumberOfDimensions==2.OR.NumberOfDimensions==3) THEN
        WRITE(14,'("    ", es25.16 )')NodeDIFFYValue(I)
      ENDIF
      IF(NumberOfDimensions==3) THEN
        WRITE(14,'("    ", es25.16 )')NodeDIFFZValue(I)
      ENDIF
      WRITE(14,'("    ", es25.16 )')NodeSTORAGEValue(I)



      IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) & 
        & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
          & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
          !source field
          IF( OUTPUT_SOURCE ) THEN
            WRITE(14,'("    ", es25.16 )')NodeSourceValue1(I)
          END IF
      END IF
    END DO
 
    WRITE(14,*) ' '
    CLOSE(14)

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

  END SUBROUTINE REACTION_DIFFUSION_IO_WRITE_NODES_CMGUI

  ! OK
  !================================================================================================================================
  !

  !> Executes element writing process.
  SUBROUTINE REACTION_DIFFUSION_IO_WRITE_ELEMENTS_CMGUI(NAME)

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
    WRITE(5,*) 'Group name: 09ryrNtatsCell'


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
        WRITE(5,*)' 2) dependent,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfVariableComponents))
      ELSE IF(I==3) THEN
        WRITE(5,*)' 3) material,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfMaterialComponents))
      ELSE IF(I==4) THEN
        WRITE(5,*)' 4) source,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfSourceComponents))
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
          SimplexOutputHelp(2)=ElementNodes(K,4)
          SimplexOutputHelp(3)=ElementNodes(K,2)
          SimplexOutputHelp(4)=ElementNodes(K,3)

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

  END SUBROUTINE REACTION_DIFFUSION_IO_WRITE_ELEMENTS_CMGUI

  ! OK
  !================================================================================================================================
  !

  !> Reorders the element node definition as needed by OpenCMISS.
  SUBROUTINE REACTION_DIFFUSION_IO_ORDER_NUMBERING(NEW,OLD,n,m,I)

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
  END SUBROUTINE REACTION_DIFFUSION_IO_ORDER_NUMBERING

  ! OK
  !================================================================================================================================
  !

  !> Writes nice information on screen.
  SUBROUTINE REACTION_DIFFUSION_IO_PRINT_ON_SCREEN

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

  END SUBROUTINE REACTION_DIFFUSION_IO_PRINT_ON_SCREEN

  ! OK
  !================================================================================================================================
  !


END MODULE REACTION_DIFFUSION_IO_ROUTINES
