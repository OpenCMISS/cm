!> \file
!> \author Chris Bradley
!> \brief This module contains all type definitions in order to avoid cyclic module references.
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

!#### Index: ne
!###  Description:
!###    Index label for a element.
!#### Index: ng
!###  Description:
!###    Index label for a gauss point.
!#### Index: ni
!###  Description:
!###    Index label for a xi direction.
!#### Index: nk
!###  Description:
!###    Index label for a derivative with respect to the global directions.
!#### Index: nn
!###  Description:
!###    Index for a local node within an element.
!#### Index: np
!###  Description:
!###    Index for a node.
!#### Index: ns
!###  Description:
!###    Index for a element parameter within an element.
!#### Index: nu
!###  Description:
!###    Index for a partial derivative.

!> This module contains all type definitions in order to avoid cyclic module references.
MODULE TYPES

  USE CMISS_PETSC_TYPES, ONLY : PETSC_ISCOLORING_TYPE,PETSC_KSP_TYPE,PETSC_MAT_TYPE,PETSC_MATFDCOLORING_TYPE,PETSC_PC_TYPE, &
    & PETSC_SNES_TYPE,PETSC_VEC_TYPE
  USE CONSTANTS
  USE KINDS
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
  USE TREES
  use linkedlist_routines
  IMPLICIT NONE

  !
  !================================================================================================================================
  !
  ! Base types
  !

  TYPE REAL_DP_PTR_TYPE
    REAL(DP), POINTER :: PTR(:)
  END TYPE REAL_DP_PTR_TYPE
  
  TYPE INTEGER_INTG_PTR_TYPE
    INTEGER(INTG), POINTER :: PTR(:)
  END TYPE INTEGER_INTG_PTR_TYPE

  TYPE INTEGER_CINT_ALLOC_TYPE
    INTEGER(C_INT), ALLOCATABLE :: ARRAY(:)
  END TYPE INTEGER_CINT_ALLOC_TYPE
  
  !
  !================================================================================================================================
  !
  ! List types
  !

  !>Buffer type to allow arrays of pointers to a list
  TYPE LIST_PTR_TYPE
    TYPE(LIST_TYPE), POINTER :: PTR !<The pointer to the list
  END TYPE LIST_PTR_TYPE

  !>Contains information on a list
  TYPE LIST_TYPE
    LOGICAL :: MUTABLE !<If true, list entries can be changed once they have been added. False by default.
    LOGICAL :: LIST_FINISHED !<Is .TRUE. if the list has finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_IN_LIST !<The number of items currently in the list
    INTEGER(INTG) :: DATA_DIMENSION !<The dimension of the data being stored
    INTEGER(INTG) :: INITIAL_SIZE !<The size of the list when it was initially created.
    INTEGER(INTG) :: SIZE !<The current size of the list.
    INTEGER(INTG) :: DATA_TYPE !<The data type of the list \see LISTS_DataType
    INTEGER(INTG) :: KEY_DIMENSION !<The key dimension number i.e., the dimension index used for indexing and sorting
    INTEGER(INTG) :: SORT_ORDER !<The ordering to be used when sorting the list \see LISTS_SortingOrder
    INTEGER(INTG) :: SORT_METHOD !<The sorting method to be used when sorting the list \see LISTS_SortingMethod
    INTEGER(INTG), ALLOCATABLE :: LIST_INTG(:) !<The integer data (dimension = 1) for integer lists. 
    INTEGER(INTG), ALLOCATABLE :: LIST_INTG2(:,:) !<The integer data (dimension > 1) for integer lists. 
    REAL(SP), ALLOCATABLE :: LIST_SP(:) !<The single precision data (dimension = 1)for single precision real lists. 
    REAL(SP), ALLOCATABLE :: LIST_SP2(:,:) !<The single precision data (dimension > 1) for single precision real lists. 
    REAL(DP), ALLOCATABLE :: LIST_DP(:) !<The double precision data (dimension = 1)for double precision real lists. 
    REAL(DP), ALLOCATABLE :: LIST_DP2(:,:) !<The double precision data (dimension > 1) for double precision real lists. 
    INTEGER(C_INT), ALLOCATABLE :: LIST_C_INT(:) !<The integer data (dimension = 1) for integer lists. 
    INTEGER(C_INT), ALLOCATABLE :: LIST_C_INT2(:,:) !<The integer data (dimension > 1) for integer lists. 
  END TYPE LIST_TYPE
    
  !
  !================================================================================================================================
  !
  ! Quadrature types
  !

  !>Contains information for a particular quadrature scheme. \see OPENCMISS::CMISSQuadratureSchemeType \todo Also evaluate the product of the basis functions at gauss points for speed???
  TYPE QUADRATURE_SCHEME_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the quadrature scheme in the list of quadrature schemes for a particular quadrature.
    TYPE(QUADRATURE_TYPE), POINTER :: QUADRATURE !<The pointer back to the quadrature for a particular quadrature scheme
    INTEGER(INTG) :: NUMBER_OF_GAUSS !<The number of gauss points for the quadrature scheme.
    REAL(DP), ALLOCATABLE :: GAUSS_POSITIONS(:,:) !<GAUSS_POSITIONS(nic,ng). The positions in the nic'th xi coordinate of Gauss point ng. Old CMISS name XIG(ni,ng,nb).
    REAL(DP), ALLOCATABLE :: GAUSS_WEIGHTS(:) !<GAUSS_WEIGHTS(ng). The weight applied to Gauss point ng. Old CMISS name WG(ng,nb).
    REAL(DP), ALLOCATABLE :: GAUSS_BASIS_FNS(:,:,:) !<GAUSS_BASIS_FNS(ns,nu,ng). The value of the basis functions evaluated at Gauss point ng for the nu'th derivative of the basis function associated with the ns'th element parameter. Old CMISS name PG(ns,nu,ng,nb)
    !Quadrature information at faces
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_FACE_GAUSS(:) !<NUMBER_OF_FACE_GAUSS(:) number of gauss points in each local face of the element
    REAL(DP), ALLOCATABLE :: FACE_GAUSS_BASIS_FNS(:,:,:,:) !<FACE_GAUSS_BASIS_FNS(ns,nu,ng,naf)
    REAL(DP), ALLOCATABLE :: FACE_GAUSS_POSITIONS(:,:,:) !<FACE_GAUSS_POSITIONS(nic,ng,naf)
    REAL(DP), ALLOCATABLE :: FACE_GAUSS_WEIGHTS(:,:) !<FACE_GAUSS_WEIGHTS(ng,naf)
  END TYPE QUADRATURE_SCHEME_TYPE

  !>A buffer type to allow for an array of pointers to a QUADRATURE_SCHEME_TYPE \see TYPES::QUADRATURE_SCHEME_TYPE
  TYPE QUADRATURE_SCHEME_PTR_TYPE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: PTR !<A pointer to the quadrature scheme
  END TYPE QUADRATURE_SCHEME_PTR_TYPE

  !>Contains information on the quadrature to be used for integrating a basis. \see OPENCMISS::CMISSQuadratureType
  TYPE QUADRATURE_TYPE
    INTEGER(INTG) :: TYPE !<The type of the quadrature \see BASIS_ROUTINES_QuadratureTypes
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The pointer back to the basis
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_GAUSS_XI(:) !<NUMBER_OF_GAUSS_XI(ni). For standard Gauss schemes the number of Gauss points to be used in the ni'th xi direction.
    INTEGER(INTG) :: GAUSS_ORDER !<For simplex Gauss schemes the order of the Quadrature scheme i.e., the order/dimension of the polynomial that can be integrated.
    TYPE(QUADRATURE_SCHEME_PTR_TYPE), ALLOCATABLE :: QUADRATURE_SCHEME_MAP(:) !<QUADRATURE_SCHEME_MAP(scheme_idx). The pointer map to the defined quadrature schemes. The size of array is given by BASIS_ROUTINES::BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES. If the quadrature scheme is not defined for the particular type then the array element is NULL. \see BASIS_ROUTINES_QuadratureSchemes.
    INTEGER(INTG) :: NUMBER_OF_SCHEMES !<The number of quadrature schemes defined for this quadrature
    TYPE(QUADRATURE_SCHEME_PTR_TYPE), POINTER :: SCHEMES(:) !<SCHEMES(scheme_idx). The array of pointers to the quadrature schemes defined for the basis. scheme_idx must be between 1 and QUADRATURE_TYPE::NUMBER_OF_SCHEMES.
    LOGICAL :: EVALUATE_FACE_GAUSS=.FALSE. !! \todo should this be here??
  END TYPE QUADRATURE_TYPE

  !
  !================================================================================================================================
  !
  ! Basis types
  !

  !> A buffer type to allow for an array of pointers to a BASIS_TYPE.
  TYPE BASIS_PTR_TYPE
    TYPE(BASIS_TYPE), POINTER :: PTR !<The pointer to the basis.
  END TYPE BASIS_PTR_TYPE

  !> Contains all information about a basis .
  TYPE BASIS_TYPE
    !\todo Add in different sub types for the different types of bases???
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the basis. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number for the basis i.e., the position indentifier for the list of bases defined.
    INTEGER(INTG) :: FAMILY_NUMBER !<The family number for the basis. A basis has a number of sub-bases attached which make a basis family. The main parent basis is the basis defined by the user and it will have a family number of 0. The sub-bases of the parent basis will be the line and face bases that make up the basis. These will have different family numbers.
    LOGICAL :: BASIS_FINISHED !<Is .TRUE. if the basis has finished being created, .FALSE. if not.
    LOGICAL :: HERMITE !<Is .TRUE. if the basis is a hermite basis, .FALSE. if not.
    INTEGER(INTG) :: TYPE !< The type of basis \see BASIS_ROUTINES_BasisTypes 
    INTEGER(INTG) :: NUMBER_OF_XI !<The number of xi directions for the basis.
    INTEGER(INTG) :: NUMBER_OF_XI_COORDINATES !<The number of xi coordinate directions for the basis. For Lagrange Hermite tensor product basis functions this is equal to the number of Xi directions. For simplex basis functions this is equal to the number of Xi directions + 1
    INTEGER(INTG), ALLOCATABLE :: INTERPOLATION_XI(:) !<INTERPOLATION_XI(ni). The interpolation specification used in the ni'th Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), ALLOCATABLE :: INTERPOLATION_TYPE(:) !<INTERPOLATION_TYPE(ni). The interpolation type in the nic'th Xi coordinate direction. Old CMISS name IBT(1,ni,nb) \see BASIS_ROUTINES_InterpolationTypes
    INTEGER(INTG), ALLOCATABLE :: INTERPOLATION_ORDER(:)!<INTERPOLATION_ORDER(ni). The interpolation order in the nic'th Xi coordinate direction. Old CMISS name IBT(2,ni,nb) \see BASIS_ROUTINES_InterpolationOrder 
    !Degenerate information
    LOGICAL :: DEGENERATE !<Is .TRUE. if the basis is a degenerate basis (i.e., has collapsed nodes), .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). The collapsed state of the ni'th direction. COLLAPSED_XI can be either XI_COLLAPSED, COLLAPSED_AT_XI0, COLLAPSED_AT_XI1 or NOT_COLLAPSED depending on whether or not the ni'th direction is collapsed, has a perpendicular Xi collapsed at the xi=0 end of the ni'th direction, has a perpendicular xi collapsed at the xi=1 of the ni'th direction or is not collapsed. NOTE: in old cmiss the value IBT(1,ni) = 5 or 6 was set for the ni that was collapsed. The perpendicular line/face was then stored in IBT(3,ni). For this the quadratic1 and quadratic2 type interpolation types are set on the perpendicular xi direction and the ni direction that is collapsed will have COLLAPSED_XI(ni) set to XI_COLLAPSED. BE CAREFUL WITH THIS WHEN TRANSLATING OLD CMISS CODE. Old CMISS name IBT(1,ni) ???? \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG) :: NUMBER_OF_COLLAPSED_XI !<The number of xi directions in the basis that are collapsed.
    LOGICAL, ALLOCATABLE :: NODE_AT_COLLAPSE(:) !<NODE_AT_COLLAPSE(nn). Is .TRUE. if the nn'th node of the basis is at a collapse, .FALSE. if not.
    !Quadrature
    TYPE(QUADRATURE_TYPE) :: QUADRATURE !<The quadrature schemes for the basis.
    INTEGER(INTG) :: NUMBER_OF_PARTIAL_DERIVATIVES !<The number of partial derivatives for the basis. Old CMISS name NUT(nbf)
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of local nodes in the basis. Old CMISS name NNT(nbf)
    !\todo
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_NODES_XIC(:) !<NUMBER_OF_NODES_XIC(nic). The number of local nodes in the nic'th coordinate in the basis. Old CMISS name IBT(2,ni,nb).
    INTEGER(INTG) :: NUMBER_OF_ELEMENT_PARAMETERS  !<The number of element parameters in the basis. Old CMISS name NST(nbf). 
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_DERIVATIVES !<The maximum number of derivatives at any node in the basis. Old CMISS name NKT(0,nbf)
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_DERIVATIVES(:) !<NUMBER_OF_DERIVATIVES(nn). The number of derivatives at the nn'th node in the basis. Old CMISS name NKT(nn,nbf).
    INTEGER(INTG), ALLOCATABLE :: NODE_POSITION_INDEX(:,:) !<NODE_POSITION_INDEX(nn,nic). The index of the node position for the nn'th local node in the nic'th coordinate. For Lagrange-Hermite tensor product basis functions: The number of coordinates equals the number of xi directions. Thus if NODE_POSITION_INDEX(nn,:)=1,2,2 then local node nn is the first node in the ni(c)=1 direction, the second node in the ni(c)=2 direction and the second node in the ni(c)=3 direction; For simplex basis functions: The number of coordinates equals the number of xi directions plus one. The index specifies the inverse distance away from the corner/end of that area coordinate. Thus if an element has quadratic interpolation the index will range from 3 (closest to the corner/end of the element that the area coordinate has the value 1) to 1 (closest to the corners/end of the element that the area coordinate has the value 0). If M is the order of the element then in NODE_POSITION_INDEX(nn,:)=1,1,M then that node is apex for the third area coordinate. In general the index values will add up to M+number of xi directions+1 (i.e., subtract one from the indicies to get the standard simplex coordinates. Old CMISS name INP(nn,ni,nb). \see TYPES::BASIS_TYPE::NODE_POSITION_INDEX_INV.
    INTEGER(INTG), ALLOCATABLE :: NODE_POSITION_INDEX_INV(:,:,:,:) !<NODE_POSITION_INDEX_INV(nnc1,nnc2,nnc3,nnc4). The inverse of the node position index for the basis. The NODE_POSITION_INDEX_INV gives the local node number for the node that has node position indices of nnc1 in the 1st ni(c) direction, nnc2 in the 2nd ni(c) direction, nnc3 in the 3rd ni(c) direction and nnc4 in the 4th ni(c) direction. NOTE: that if the basis has less than 4 ni(c) direction the position index is 1. Old CMISS name NNB(inp1,inp2,inp3,nbf). \see TYPES::BASIS_TYPE::NODE_POSITION_INDEX.
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_ORDER_INDEX(:,:,:) !<DERIVATIVE_ORDER_INDEX(nk,nn,0:ni). The index of the derivative order for the nk'th derivative of the nn'th node in the ni'th direction of the basis. The derivative index is NO_PART_DERIV for zeroth order, FIRST_PART_DERIV for the first order and SECOND_PART_DERIV for the second order derivative. Thus a DERIVATIVE_ORDER_INDEX(nk,nn,1..) of {NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV} indicates that the nk'th derivative of the nn'th node of the basis is the first derivative with respect to the s2 direction. Old CMISS name IDO(nk,nn,1:ni,nbf). \see TYPES::BASIS_TYPE::DERIVATIVE_ORDER_INDEX_INV,CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_ORDER_INDEX_INV(:,:,:,:) !<DERIVATIVE_ORDER_INDEX_INV(nu1,nu2,nu3,nn). The inverse of the derivative order index for the nn'th local node of the basis. DERIVATIVE_ORDER_INDEX_INV gives the derivative number for the nu1 partial derivative in the 1st xi direction, the nu2 partial derivative in the 2nd xi direction and the nu3 partial derivative in the 3rd xi direction. NOTE: local node nn does not carry any derivatives of the requested partial derivative type then DERIVATIVE_ORDER_INDEX_INV will return 0. If the basis has less than 3 xi directions then the nu index is 1. \see TYPES::BASIS_TYPE::DERIVATIVE_ORDER_INDEX
    INTEGER(INTG), ALLOCATABLE :: PARTIAL_DERIVATIVE_INDEX(:,:) !<PARTIAL_DERIVATIVE_INDEX(nk,nn). Gives the partial derivative number (nu) of the nk'th derivative of the nn'th local node for the basis. Old CMISS name IDO(nk,nn,0,nbf).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_PARAMETER_INDEX(:,:) !<ELEMENT_PARAMETER_INDEX(nk,nn). Gives the element parameter number (ns) of the nk'th derivative of the nn'th local node for the basis. Old CMISS name NSB(nk,nn,nbf).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_PARAMETER_INDEX_INV(:,:) !<ELEMENT_PARAMETER_INDEX_INV(1..2,ns). Gives the inverse fo the element parameter index. ELEMENT_PARAMETER_INDEX_INV(1,ns) gives the local node number corresponding to the ns'th element parameter. ELEMENT_PARAMETER_INDEX_INV(2,ns) gives the local derivative number corresponding to the ns'th element parameter.
    !Line information
    INTEGER(INTG) :: NUMBER_OF_LOCAL_LINES !<The number of local lines in the basis.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_LINE_XI_DIRECTION(:) !<LOCAL_LINE_XI_DIRECTION(nae). The Xi direction of the nae'th local line for the basis.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_NODES_IN_LOCAL_LINE(:) !<NUMBER_OF_NODES_IN_LOCAL_LINE(nae). The the number of nodes in the nae'th local line for the basis. Old CMISS name NNL(0,nae,nb).
    INTEGER(INTG), ALLOCATABLE :: NODE_NUMBERS_IN_LOCAL_LINE(:,:) !<NODE_NUMBERS_IN_LOCAL_LINE(nnl,nae). The local node numbers (nn) for the nnl'th line node in the nae'th local line for the basis. Old CMISS name NNL(1..,nae,nb).
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_NUMBERS_IN_LOCAL_LINE(:,:) !<DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nnl,nae). The derivative numbers (nk) for the nnl'th line node in the nae'th local line for the basis.
    !Face information
    INTEGER(INTG) :: NUMBER_OF_LOCAL_FACES !<The number of local faces in the basis.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_FACE_XI_DIRECTION(:) !<LOCAL_FACE_XI_DIRECTION(nae). The Xi direction of the nae'th local face for the basis.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_NODES_IN_LOCAL_FACE(:) !<NUMBER_OF_NODES_IN_LOCAL_FACE(nae). The the number of nodes in the nae'th local face for the basis. Old CMISS name NNL(0,nae,nb).
    INTEGER(INTG), ALLOCATABLE :: NODE_NUMBERS_IN_LOCAL_FACE(:,:) !<NODE_NUMBERS_IN_LOCAL_FACE(nnl,nae). The local element node numbers (nn) for the nnl'th face node in the nae'th local face for the basis. Old CMISS name NNL(1..,nae,nb).
    !\todo Remove local face node derivative below since all nodes of the face would have the same number of derivatives (unless a node is collapsed?)
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_NUMBERS_IN_LOCAL_FACE(:,:,:) !<DERIVATIVES_NUMBERS_IN_LOCAL_FACE(nnk,nnl,nae). The element derivative numbers (nk) for the nnk'th face derivative's nnl'th face node in the nae'th local face for the basis.
    !\todo What is the difference between LOCAL_XI_NORMAL and LOCAL_FACE_XI_DIRECTION ? They're the same
    INTEGER(INTG), ALLOCATABLE :: LOCAL_XI_NORMAL(:) !<LOCAL_XI_NORMAL(nae). The Xi direction that is normal to either the nae'th local line for bases with 2 xi directions or the nae'th local face for bases with 3 xi directions. For bases with 1 xi direction the array is not allocated. Note: Normals are always outward.
    !Sub-basis information
    TYPE(BASIS_PTR_TYPE), POINTER :: LINE_BASES(:) !<LINE_BASES(nae). The pointer to the basis for the nae'th line for the basis.
    TYPE(BASIS_PTR_TYPE), POINTER :: FACE_BASES(:) !<FACE_BASES(naf). The pointer to the basis for the naf'th face for the basis.
    INTEGER(INTG) :: NUMBER_OF_SUB_BASES !<The number of sub-bases (lines, faces) for the basis.
    TYPE(BASIS_PTR_TYPE), POINTER :: SUB_BASES(:) !<SUB_BASES(sbn). The pointer to the sbn'th sub-basis for the basis.
    TYPE(BASIS_TYPE), POINTER :: PARENT_BASIS !<The pointer to the parent basis for the basis. NOTE: that if the basis is not a sub-basis of another basis this pointer will be NULL. 
  END TYPE BASIS_TYPE

  !>Contains information on the defined basis functions
  TYPE BASIS_FUNCTIONS_TYPE
   INTEGER(INTG) :: NUMBER_BASIS_FUNCTIONS !<The number of basis functions defined
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:) !<The array of pointers to the defined basis functions
  END TYPE BASIS_FUNCTIONS_TYPE
  
  !
  !================================================================================================================================
  !
  ! Coordinate system types
  !

  !>Contains information on a coordinate system. \todo Have a list of coordinate systems and have a pointer in the coordinate_system_type back to the regions that use them.
  TYPE, BIND(C) :: COORDINATE_SYSTEM_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the coordinate. The user number must be unique.
    LOGICAL :: COORDINATE_SYSTEM_FINISHED !<Is .TRUE. if the coordinate system has finished being created, .FALSE. if not.
    INTEGER(INTG) :: TYPE !<The type of coordinate system. Old CMISS name ITYP10(nr). \see COORINDATE_ROUTINES_CoordinateSystemTypes
    INTEGER(INTG) :: RADIAL_INTERPOLATION_TYPE !<The type of radial interpolation type for non-rectangular cartesian systems. Old CMISS name JTYP10(nr). \see COORDINATE_ROUTINES_RadialInterpolations
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS !<The number of dimensions for the coordinate system. Old CMISS name NJT.
    REAL(DP) :: FOCUS !<The focus of the coordinate system for a prolate-spheriodal coordinate system.
    REAL(DP) :: ORIGIN(3) !<ORIGIN(nj). The nj'th component of the origin of the coordinate system wrt the global coordinate system. NOTE: maybe this should be wrt to the parent regions coordinate system - this would then go into the REGION type.
    REAL(DP) :: ORIENTATION(3,3) !<ORIENTATION(nj,mj). he orientation matrix for the orientation of the coordinate system wrt to the global coordinate system. NOTE: maybe this should be wrt to the parent regions coordinate system - this would then go into the REGION type.
  END TYPE COORDINATE_SYSTEM_TYPE

  !>A buffer type to allow for an array of pointers to a COORDINATE_SYSTEM_TYPE.
  TYPE COORDINATE_SYSTEM_PTR_TYPE
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: PTR !<A pointer to the coordinate system
  END TYPE COORDINATE_SYSTEM_PTR_TYPE
  
  !
  !================================================================================================================================
  !
  ! Data point types
  !

  !>Contains information about a data point.
  TYPE DATA_POINT_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of data point. 
    INTEGER(INTG) :: USER_NUMBER !<The user defined number of data point. 
    TYPE(VARYING_STRING) :: LABEL !<A string label for the data point.
    REAL(DP), ALLOCATABLE :: VALUES(:) !Values of the data point specifying the spatial position in the region, has the size of region dimension the data point belongs to.
    REAL(DP), ALLOCATABLE :: WEIGHTS(:) !<Weights of the data point, has the size of region dimension the data point belongs to.
    REAL(DP) :: PROJECTION_DISTANCE !<The distances between the data point and the projection. Assigned only if DATA_POINTS_PROJECTED is .TRUE.
    INTEGER(INTG) :: PROJECTION_ELEMENT_NUMBER !<The element of the mesh the data point projects onto. Assigned only if DATA_POINTS_PROJECTED is .TRUE.
    INTEGER(INTG) :: PROJECTION_ELEMENT_FACE_NUMBER !<The element face of the mesh the data point projects onto. Assigned only if DATA_POINTS_PROJECTED is .TRUE. and DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE is chosen    
    INTEGER(INTG) :: PROJECTION_ELEMENT_LINE_NUMBER !<The element line of the mesh the data point projects onto. Assigned only if DATA_POINTS_PROJECTED is .TRUE. and DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE is chosen
    INTEGER(INTG) :: PROJECTION_EXIT_TAG !<The exit tage of the data projection. \See DATA_PROJECTION_ROUTINES. Assigned only if DATA_POINTS_PROJECTED is .TRUE. 
    REAL(DP), ALLOCATABLE :: PROJECTION_XI(:) !<The xi coordinate of the projection. Assigned only if DATA_POINTS_PROJECTED is .TRUE.
  END TYPE DATA_POINT_TYPE

  !>Contains information on the data points defined on a region. \see OPENCMISS::CMISSDataPointsType
  TYPE DATA_POINTS_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the data points. If the data points are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the data points. If the data points are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    LOGICAL :: DATA_POINTS_FINISHED !<Is .TRUE. if the data points have finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_DATA_POINTS !<The number of data points defined on the region.
    TYPE(DATA_POINT_TYPE), ALLOCATABLE :: DATA_POINTS(:) !<DATA_POINTS(data_points_idx). The data point information for the data_points_idx'th global data point.
    TYPE(TREE_TYPE), POINTER :: DATA_POINTS_TREE !<The tree for user to global data point mapping.
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for the data points.
    LOGICAL :: DATA_POINTS_PROJECTED !<Is .TRUE. if the data points have been projected, .FALSE. if not.
  END TYPE DATA_POINTS_TYPE
  
  !
  !================================================================================================================================
  !
  ! Data projection types
  !

  TYPE DATA_PROJECTION_TYPE
    LOGICAL :: DATA_PROJECTION_FINISHED !<Is .TRUE. if the data projection has finished being created, .FALSE. if not.
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<The pointer to the data points for this data projection.
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<The pointer to the geometric field for this data projection.
    INTEGER(INTG) :: COORDINATE_SYSTEM_DIMENSIONS !<The coordinate system dimension of this data projection.
    REAL(DP) :: MAXIMUM_ITERATION_UPDATE !<The maximum xi update allowed at each newton iteration, analogous to maximum trust region size in the trust region model approach.
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ITERATIONS !<The maximum number of iterations
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_ELEMENTS !<The number of closest elements to perform full projection on. The algorithm first find the distance of the data point to each elements base on starting xi, full projection is only performed on the first few elements sorted by the distance
    INTEGER(INTG) :: NUMBER_OF_XI !<The number of xi of the mesh, ie. the mesh dimension
    INTEGER(INTG) :: PROJECTION_TYPE !<type of projection to perform. \See DATA_PROJECTION_ROUTINES     
    REAL(DP), ALLOCATABLE :: STARTING_XI(:) !<The starting value of the element xi
    REAL(DP) :: ABSOLUTE_TOLERANCE !<The absolute tolerance of the iteration update
    REAL(DP) :: RELATIVE_TOLERANCE !<The relative tolerance of the iteration update
  END TYPE DATA_PROJECTION_TYPE

  !
  !================================================================================================================================
  !
  ! Node types
  !

  !>Contains information about a node.
  TYPE NODE_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the node.
    INTEGER(INTG) :: USER_NUMBER !<The user defined number of the node.
    TYPE(VARYING_STRING) :: LABEL !<A string label for the node
  END TYPE NODE_TYPE

  !>Contains information on the nodes defined on a region. \see OPENCMISS::CMISSNodesType
  TYPE NODES_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the nodes. If the nodes are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the nodes. If the nodes are in a region rather than an interface then this pointer will be NULL and the region pointer should be used.
    LOGICAL :: NODES_FINISHED !<Is .TRUE. if the nodes have finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of nodes defined on the region.
    TYPE(NODE_TYPE), ALLOCATABLE :: NODES(:) !<NODES(nodes_idx). The nodal information for the nodes_idx'th global node.
    TYPE(TREE_TYPE), POINTER :: NODES_TREE !<The tree for user to global node mapping
  END TYPE NODES_TYPE

  !
  !================================================================================================================================
  !
  ! Mesh types
  !

  !>Contains information on the dofs for a mesh.
  TYPE MESH_DOFS_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh.
    INTEGER(INTG) :: NUMBER_OF_DOFS !<The number of dofs in the mesh.
  END TYPE MESH_DOFS_TYPE

  !>Contains information on the mesh adjacent elements for a xi coordinate
  TYPE MESH_ADJACENT_ELEMENT_TYPE
    INTEGER(INTG) :: NUMBER_OF_ADJACENT_ELEMENTS !<The number of adjacent elements for the xi coordinate
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:) !<The global element numbers of the elements adjacent to this element for the xi coordinate
  END TYPE MESH_ADJACENT_ELEMENT_TYPE
  
  !>Contains the information for an element in a mesh.
  TYPE MESH_ELEMENT_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global element number in the mesh.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the element.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function for the element.
    INTEGER(INTG), ALLOCATABLE :: MESH_ELEMENT_NODES(:) !<MESH_ELEMENT_NODES(nn). The mesh node number in the mesh of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:) !<GLOBAL_ELEMENT_NODES(nn). The global node number in the mesh of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: USER_ELEMENT_NODE_VERSIONS(:,:) !< USER_ELEMENT_NODE_VERSIONS(derivative_idx,element_node_idx).  The version number for the nn'th user node's nk'th derivative. Size of array dependent on the maximum number of derivatives for the basis of the specified element.
    INTEGER(INTG), ALLOCATABLE :: USER_ELEMENT_NODES(:) !<USER_ELEMENT_NODES(nn). The user node number in the mesh of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    TYPE(MESH_ADJACENT_ELEMENT_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS(:) !<ADJACENT_ELEMENTS(-nic:nic). The adjacent elements information in the nic'th xi coordinate direction. Note that -nic gives the adjacent element before the element in the nic'th direction and +nic gives the adjacent element after the element in the nic'th direction. The ni=0 index will give the information on the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    !INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ADJACENT_ELEMENTS(:) !<NUMBER_OF_ADJACENT_ELEMENTS(-ni:ni). The number of elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent element before the element in the ni'th direction and +ni gives the adjacent element after the element in the ni'th direction. The ni=0 index should be 1 for the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    !INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:,:) !<ADJACENT_ELEMENTS(nei,-ni:ni). The local element numbers of the elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent elements before the element in the ni'th direction and +ni gives the adjacent elements after the element in the ni'th direction. The ni=0 index should give the current element number. Old CMISS name NXI(-ni:ni,0:nei,ne)
    LOGICAL :: BOUNDARY_ELEMENT !<Is .TRUE. if the mesh element is on the boundary of a mesh, .FALSE. if not.
  END TYPE MESH_ELEMENT_TYPE

  !>Contains the information for the elements of a mesh.
  TYPE MESH_ELEMENTS_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<The pointer to the mesh for the elements information.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !< The number of elements in the mesh.
    LOGICAL :: ELEMENTS_FINISHED !<Is .TRUE. if the mesh elements have finished being created, .FALSE. if not.
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENTS(:) !<ELEMENTS(ne). The pointer to the array of information for the elements of this mesh. ELEMENTS(ne) contains the information for the ne'th global element of the mesh. \todo Should this be allocatable.
    TYPE(TREE_TYPE), POINTER :: ELEMENTS_TREE !<A tree mapping the mesh global element number to the mesh user element number.
  END TYPE MESH_ELEMENTS_TYPE

  !>Contains the information for a node derivative of a mesh.
  TYPE MESH_NODE_DERIVATIVE_TYPE
    INTEGER(INTG) :: NUMBER_OF_VERSIONS !The number of global versions at the node for the mesh.
    INTEGER(INTG), ALLOCATABLE :: USER_VERSION_NUMBERS(:) !The user version index of the nk'th global derivative for the node.
    INTEGER(INTG), ALLOCATABLE :: DOF_INDEX(:) !The global dof version index (nv) in the domain of the nk'th global derivative for the node.
    INTEGER(INTG) :: GLOBAL_DERIVATIVE_INDEX !The global derivative index of the nk'th global derivative for the node.
    INTEGER(INTG) :: PARTIAL_DERIVATIVE_INDEX !The partial derivative index (nu) of the nk'th global derivative for the node. Old CMISS name NUNK(nk,nj,np).
  END TYPE MESH_NODE_DERIVATIVE_TYPE

  !>Contains the topology information for a global node of a mesh.
  TYPE MESH_NODE_TYPE
    INTEGER(INTG) :: MESH_NUMBER !<The mesh node number in the mesh.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global node number in the mesh.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the node.
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES !<The number of global derivatives at the node for the mesh. Old CMISS name NKT(nj,np).
    TYPE(MESH_NODE_DERIVATIVE_TYPE), ALLOCATABLE :: DERIVATIVES(:) !<DERIVATIVES(derivative_idx)
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements surrounding the node in the mesh. Old CMISS name NENP(np,0,0:nr).
    INTEGER(INTG), POINTER :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nep). The global element number of the nep'th element that is surrounding the node. Old CMISS name NENP(np,nep,0:nr). \todo Change this to allocatable.
    LOGICAL :: BOUNDARY_NODE !<Is .TRUE. if the mesh node is on the boundary of the mesh, .FALSE. if not.
  END TYPE MESH_NODE_TYPE

  !>Contains the information for the nodes of a mesh.
  TYPE MESH_NODES_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<The pointer to the mesh for this nodes information.
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of nodes in the mesh.
    TYPE(MESH_NODE_TYPE), POINTER :: NODES(:) !<NODES(np). The pointer to the array of topology information for the nodes of the mesh. NODES(np) contains the topological information for the np'th global node of the mesh. \todo Should this be allocatable???
    TYPE(TREE_TYPE), POINTER :: NODES_TREE !<A tree mapping the mesh global number to the region nodes global number.
  END TYPE MESH_NODES_TYPE

  !>Contains information on the (global) topology of a mesh.
  TYPE MESH_TOPOLOGY_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<Pointer to the parent mesh.
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The mesh component number for this mesh topology.
    TYPE(MESH_NODES_TYPE), POINTER :: NODES !<Pointer to the nodes within the mesh topology.
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<Pointer to the elements within the mesh topology.
    TYPE(MESH_DOFS_TYPE), POINTER :: DOFS !<Pointer to the dofs within the mesh topology.
  END TYPE MESH_TOPOLOGY_TYPE

  !>A buffer type to allow for an array of pointers to a MESH_TOPOLOGY_TYPE.
  TYPE MESH_TOPOLOGY_PTR_TYPE
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: PTR !<The pointer to the mesh topology.
  END TYPE MESH_TOPOLOGY_PTR_TYPE

  !>Embedded mesh types
  TYPE EMBEDDING_XI_TYPE
    INTEGER(INTG) :: NUMBER_OF_NODES !<Number of nodes embedded in this element
    INTEGER(INTG), ALLOCATABLE :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in this element
    REAL(DP), ALLOCATABLE :: XI_COORDS(:,:) !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node this element
  END TYPE EMBEDDING_XI_TYPE

  TYPE EMBEDDING_GAUSSPOINT_TYPE
    INTEGER(INTG) :: ELEMENT_NUMBER !<Element number in child mesh
    REAL(DP), ALLOCATABLE :: CHILD_XI_COORD(:) !<Xi coord in this element
    REAL(DP), ALLOCATABLE :: PARENT_XI_COORD(:) !<Xi coordinates in parent element, not really needed but can be useful
  END TYPE EMBEDDING_GAUSSPOINT_TYPE

  TYPE MESH_EMBEDDING_TYPE
    TYPE(MESH_TYPE), POINTER :: CHILD_MESH, PARENT_MESH !<The parent and child mesh
    TYPE(EMBEDDING_XI_TYPE), ALLOCATABLE :: CHILD_NODE_XI_POSITION(:)            !<Location of embedded nodes for each element in the parent mesh
    TYPE(EMBEDDING_GAUSSPOINT_TYPE), ALLOCATABLE :: GAUSS_POINT_XI_POSITION(:,:) !<GAUSS_POINT_XI_POSITION(gauss_idx,element_idx) Location of the gauss_idx'th Gauss point of the element_idx'th parent mesh in the child mesh
  END TYPE MESH_EMBEDDING_TYPE


  !>Contains information on a mesh defined on a region. \see OPENCMISS::CMISSMeshType
  TYPE MESH_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user number of the mesh. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global number for the mesh.
    LOGICAL :: MESH_FINISHED !<Is .TRUE. if the mesh has finished being created, .FALSE. if not.
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes for this mesh.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the mesh. If the mesh is in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the mesh. If the mesh is in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh generate this mesh.
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS !<The number of dimensions (Xi directions) for this mesh.
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !<The number of mesh components in this mesh.
    LOGICAL :: MESH_EMBEDDED !<Is .TRUE. if the mesh is embedded in another mesh, .FALSE. if not.
    TYPE(MESH_TYPE), POINTER :: EMBEDDING_MESH !<If this mesh is embedded the pointer to the mesh that this mesh is embedded in. IF the mesh is not embedded the pointer is NULL.
    INTEGER(INTG) :: NUMBER_OF_EMBEDDED_MESHES !<The number of meshes that are embedded in this mesh.
    TYPE(MESH_PTR_TYPE), POINTER :: EMBEDDED_MESHES(:) !<EMBEDDED_MESHES(mesh_idx). A pointer to the mesh_idx'th mesh that is embedded in this mesh.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS
    TYPE(MESH_TOPOLOGY_PTR_TYPE), POINTER :: TOPOLOGY(:) !<TOPOLOGY(mesh_component_idx). A pointer to the topology mesh_component_idx'th mesh component.
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS !<A pointer to the decompositions for this mesh.
    LOGICAL :: SURROUNDING_ELEMENTS_CALCULATE !<Boolean flag to determine whether surrounding elements should be calculated.
  END TYPE MESH_TYPE

  !>A buffer type to allow for an array of pointers to a MESH_TYPE.
  TYPE MESH_PTR_TYPE
    TYPE(MESH_TYPE), POINTER :: PTR !<The pointer to the mesh. 
  END TYPE MESH_PTR_TYPE

  !>Contains information on the meshes defined on a region.
  TYPE MESHES_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containg the meshes. If the meshes are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containg the meshes. If the meshes are in a region rather than an interface then this pointer will be NULL and the region pointer should be used.
    INTEGER(INTG) :: NUMBER_OF_MESHES !<The number of meshes defined on the region.
    TYPE(MESH_PTR_TYPE), POINTER :: MESHES(:) !<MESHES(meshes_idx). The array of pointers to the meshes.
  END TYPE MESHES_TYPE

  !
  !================================================================================================================================
  !
  ! Generated Mesh types
  !

  !>Contains information on a generated regular mesh
  TYPE GENERATED_MESH_REGULAR_TYPE
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: BASES(:) !<The pointers to the bases used in the regular mesh.
    INTEGER(INTG) :: COORDINATE_DIMENSION !<The number of coordinates for the regular mesh.
    INTEGER(INTG) :: MESH_DIMENSION !<The dimension/number of Xi directions of the regular mesh.
    REAL(DP), ALLOCATABLE :: ORIGIN(:) !<ORIGIN(coordinate_idx). The position of the origin (first) corner of the regular mesh
    REAL(DP), ALLOCATABLE :: MAXIMUM_EXTENT(:) !<MAXIMUM_EXTENT(coordinate_idx). The extent/size in each nj'th direction of the regular mesh.
    REAL(DP), ALLOCATABLE :: BASE_VECTORS(:,:) !<MESH_BASE_VECTORS(coordinate_idx,xi_idx). The base vectors indicating the geometric direction of the xi_idx mesh coordinate.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ELEMENTS_XI(:) !<NUMBER_OF_ELEMENTS_XI(xi_idx). The number of elements in the xi_idx'th Xi direction for the mesh.
  END TYPE GENERATED_MESH_REGULAR_TYPE

  !>Contains information of a generated cylinder mesh
  !>Allows only a 3D cylinder mesh with xi directions (r,theta,z)
  TYPE GENERATED_MESH_CYLINDER_TYPE
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh.
    REAL(DP), ALLOCATABLE :: ORIGIN(:) !<ORIGIN(nj). The position of the origin (centre) of lower face of cylinder mesh.
    REAL(DP), ALLOCATABLE :: CYLINDER_EXTENT(:) !<CYLINDER_EXTENT(nj). The size of inner & outer radii and height of cylinder.
    INTEGER(INTG) :: MESH_DIMENSION !<The dimension/number of Xi directions of the cylinder mesh.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ELEMENTS_XI(:) !<NUMBER_OF_ELEMENTS(ni). The number of elements in radial, circumferential and axial directions
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: BASES(:) !<The pointers to the bases used in the regular mesh.
    LOGICAL :: APPEND_LINEAR_COMPONENT=.FALSE. !<True when two mesh components are needed
 END TYPE GENERATED_MESH_CYLINDER_TYPE
 
  !>Contains information of a generated ellipsoid mesh
  !>Allows only a 3D ellipsoid mesh
 TYPE GENERATED_MESH_ELLIPSOID_TYPE
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh.
    REAL(DP), ALLOCATABLE :: ORIGIN(:) !<ORIGIN(nj). The position of the origin (centre) of lower face of ellipsoid mesh.
    REAL(DP), ALLOCATABLE :: ELLIPSOID_EXTENT(:) !<ELLIPSOID_EXTENT(nj). The size of long axis, short axis, wall thickness and cut off angle of ellipsoid.
    INTEGER(INTG) :: MESH_DIMENSION !<The dimension/number of Xi directions of the ellipsoid mesh.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ELEMENTS_XI(:) !<NUMBER_OF_ELEMENTS(ni). The number of elements in circumferential, longitudinal and transmural directions
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: BASES(:) !<The pointers to the bases used in the ellipsoid mesh
    LOGICAL :: APPEND_LINEAR_COMPONENT=.FALSE. !<True when two mesh components are needed 
END TYPE GENERATED_MESH_ELLIPSOID_TYPE

  !>Contains information on a generated mesh. \see OPENCMISS::CMISSGeneratedMeshType
  TYPE GENERATED_MESH_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user number of the generated mesh. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global number for the generated mesh.
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<Is .TRUE. if the generated mesh has finished being created, .FALSE. if not.
    LOGICAL :: GENERATED_MESH_FINISHED !<Is .TRUE. if the generated mesh has finished being created, .FALSE. if not.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the generated mesh. If the generated mesh is in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the generated mesh. If the generated mesh is in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: GENERATED_TYPE !<The type of the generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH !<A pointer to the regular generated mesh information if the generated mesh is a regular mesh, NULL if not.
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH !<A pointer to the cylinder generate mesh information if the generated mesh is a cylinder mesh, NULL if not.
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH !<A pointer to the ellipsoid generate mesh information if the generated mesh is a ellipsoid mesh, NULL if not.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh that has been generated.
  END TYPE GENERATED_MESH_TYPE
  
  !>A buffer type to allow for an array of pointers to a GENERATED_MESH_TYPE.
  TYPE GENERATED_MESH_PTR_TYPE
    TYPE(GENERATED_MESH_TYPE), POINTER :: PTR !<The pointer to the generated mesh.
  END TYPE GENERATED_MESH_PTR_TYPE
       
  !>Contains information on the generated meshes defined on a region.
  TYPE GENERATED_MESHES_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the generated meshes. If the generated meshes are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the generated meshes. If the generated meshes are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: NUMBER_OF_GENERATED_MESHES !<The number of generated meshes defined.
    TYPE(GENERATED_MESH_PTR_TYPE), POINTER :: GENERATED_MESHES(:) !<The array of pointers to the generated meshes.
   END TYPE GENERATED_MESHES_TYPE
  
  !
  !================================================================================================================================
  !
  ! Domain types
  !
  
  !>Contains information on the degrees-of-freedom (dofs) for a domain.
  TYPE DOMAIN_DOFS_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain.
    INTEGER(INTG) :: NUMBER_OF_DOFS !<The number of degrees-of-freedom (excluding ghost dofs) in the domain.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_DOFS !<The total number of degrees-of-freedom (including ghost dofs) in the domain.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_DOFS !<The number of global degrees-of-freedom in the domains.
    INTEGER(INTG), ALLOCATABLE :: DOF_INDEX(:,:) !<DOF_INDEX(i,ny). The index for the ny'th degree-of-freedom. When i=1 DOF_INDEX will give the global derivative version number (version_idx) associated with the dof. When i=2 DOF_INDEX will give the global derivative number (derivative_idx) associated with the dof. When i=3 DOF_INDEX will give the local node number (node_idx) associated with the dof.
  END TYPE DOMAIN_DOFS_TYPE

  !>Contains the information for a line in a domain.
  TYPE DOMAIN_LINE_TYPE
    INTEGER(INTG) :: NUMBER !<The line number in the domain.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function for the line.
    INTEGER(INTG), ALLOCATABLE :: NODES_IN_LINE(:) !<NODES_IN_LINE(nn). The local node number in the domain of the nn'th local node in the line. Old CMISS name NPL(2..5,nj,nl).
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVES_IN_LINE(:,:,:) !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx). When i=1 DERIVATIVES_IN_LINE will give the global derivative number of the local_derivative_idx'th local derivative for the local_node_idx'th local node in the line, When i=2 DERIVATIVES_IN_LINE will give the global derivative version number of the local_derivative_idx'th local derivative for the local_node_idx'th local node in the line. Old CMISS name NPL(4..5,nj,nl).
    LOGICAL :: BOUNDARY_LINE !<Is .TRUE. if the line is on the boundary of the mesh for the domain, .FALSE. if not.
    INTEGER(INTG) :: ELEMENT_NUMBER !<The element number of the element on which the line is on
  END TYPE DOMAIN_LINE_TYPE

  !>A buffer type to allow for an array of pointers to a DOMAIN_LINE_TYPE
  TYPE DOMAIN_LINE_PTR_TYPE
    TYPE(DOMAIN_LINE_TYPE), POINTER :: PTR !<A pointer to the domain line.
  END TYPE DOMAIN_LINE_PTR_TYPE

  !>Contains the topology information for the lines of a domain.
  TYPE DOMAIN_LINES_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this lines topology information.
    INTEGER(INTG) :: NUMBER_OF_LINES !<The number of lines in this domain topology.
    TYPE(DOMAIN_LINE_TYPE), ALLOCATABLE :: LINES(:) !<LINES(nl). The pointer to the array of topology information for the lines of this domain. LINES(nl) contains the topological information for the nl'th local line of the domain.
  END TYPE DOMAIN_LINES_TYPE

  !>Contains the information for a face in a domain.
  TYPE DOMAIN_FACE_TYPE
    INTEGER(INTG) :: NUMBER !<The face number in the domain.
    INTEGER(INTG) :: XI_DIRECTION1 !<The first xi direction of the face. \todo move this to the decomposition face type
    INTEGER(INTG) :: XI_DIRECTION2 !<The second xi direction of the face. \todo move this to the decomposition face type
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function for the face.
    INTEGER(INTG), ALLOCATABLE :: NODES_IN_FACE(:) !<NODES_IN_FACE(nn). The local node number in the domain of the nn'th local node in the face. Old CMISS name NPNF(nn,nbf).
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVES_IN_FACE(:,:,:) !<DERIVATIVES_IN_FACE(i,local_derivative_idx,local_node_idx). When i=1 DERIVATIVES_IN_FACE will give the global derivative number of the local derivative_idx'th local derivative for the local_node_idx'th local node in the face, When i=2 DERIVATIVES_IN_FACE will give the global derivative version number of the local derivative_idx'th local derivative for the local_node_idx'th local node in the face.
    LOGICAL :: BOUNDARY_FACE !<Is .TRUE. if the face is on the boundary of the mesh for the domain, .FALSE. if not.
    INTEGER(INTG) :: ELEMENT_NUMBER !<The element number of the element on which the face is on
  END TYPE DOMAIN_FACE_TYPE

  !>A buffer type to allow for an array of pointers to a DOMAIN_FACE_TYPE.
  TYPE DOMAIN_FACE_PTR_TYPE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: PTR !<The pointer to the domain face.
  END TYPE DOMAIN_FACE_PTR_TYPE

  !>Contains the topology information for the faces of a domain.
  TYPE DOMAIN_FACES_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this faces topology information.
    INTEGER(INTG) :: NUMBER_OF_FACES !<The number of faces in this domain topology.
    TYPE(DOMAIN_FACE_TYPE), ALLOCATABLE :: FACES(:) !<FACES(nf). The pointer to the array of topology information for the faces of this domain. FACES(nf) contains the topological information for the nf'th local face of the domain.
  END TYPE DOMAIN_FACES_TYPE

  !>Contains the information for an element in a domain.
  TYPE DOMAIN_ELEMENT_TYPE
    INTEGER(INTG) :: NUMBER !<The local element number in the domain.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function for the element.
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:) !<ELEMENT_NODES(element_node_idx). The local node number in the domain of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_DERIVATIVES(:,:,:) !<ELEMENT_DERIVATIVES(i,local_derivative_idx,local_element_node_idx).  When i=1 ELEMENT_DERIVATIVES will give the global derivative number of the local_derivative_idx'th local derivative for the local_element_node_idx'th local node in the element. When i=2 ELEMENT_DERIVATIVES will give the version number of the local_derivative_idx'th local derivative for the local_element_node_idx'th local node in the element. Old CMISS name NKJE(nk,nn,nj,ne).
  END TYPE DOMAIN_ELEMENT_TYPE
  
  !>Contains the topology information for the elements of a domain.
  TYPE DOMAIN_ELEMENTS_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this elements topology information.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !<The number of elements (excluding ghost elements) in this domain topology.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ELEMENTS !<The total number of elements (including ghost elements) in this domain topology.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ELEMENTS !<The number of global elements in this domain topology.
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: ELEMENTS(:) !<ELEMENTS(ne). The pointer to the array of topology information for the elements of this domain. ELEMENTS(ne) contains the topological information for the ne'th local elements of the domain. \todo Change this to allocatable???
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS !<The maximum number of element parameters (ns) for all the elements in the domain.
  END TYPE DOMAIN_ELEMENTS_TYPE

  !>Contains the topology information for a local node derivative of a domain.
  TYPE DOMAIN_NODE_DERIVATIVE_TYPE
    INTEGER(INTG) :: NUMBER_OF_VERSIONS !The number of global versions at the node for the mesh.
    INTEGER(INTG), ALLOCATABLE :: USER_VERSION_NUMBERS(:) !The user version index of the nk'th global derivative for the node.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_VERSION_NUMBERS(:) !The local version index of the local nk'th global derivative for the node in the domain.
    INTEGER(INTG), ALLOCATABLE :: DOF_INDEX(:) !The local dof derivative version index in the domain of the nk'th global derivative for the node.
    INTEGER(INTG) :: GLOBAL_DERIVATIVE_INDEX !The global derivative index of the nk'th global derivative for the node.
    INTEGER(INTG) :: PARTIAL_DERIVATIVE_INDEX !The partial derivative index (nu) of the nk'th global derivative for the node. Old CMISS name NUNK(nk,nj,np).
  END TYPE DOMAIN_NODE_DERIVATIVE_TYPE

  !>Contains the topology information for a local node of a domain.
  TYPE DOMAIN_NODE_TYPE
    INTEGER(INTG) :: LOCAL_NUMBER !<The local node number in the domain.
    INTEGER(INTG) :: MESH_NUMBER !<The corresponding global node number in the mesh of the local node number in the domain i.e., the mesh node number.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global number for the node i.e., the node number in the list of nodes for the region.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the node.
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES !<The number of global derivatives at the node for the domain. Old CMISS name NKT(nj,np)
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE), ALLOCATABLE :: DERIVATIVES(:) !<DERIVATIVES(derivative_idx)
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements surrounding the node in the domain. Old CMISS name NENP(np,0,0:nr).
    INTEGER(INTG), POINTER :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nep). The local element number of the nep'th element that is surrounding the node. Old CMISS name NENP(np,nep,0:nr). \todo Change this to allocatable.
    INTEGER(INTG) :: NUMBER_OF_NODE_LINES !<The number of lines surrounding the node in the domain.
    INTEGER(INTG), ALLOCATABLE :: NODE_LINES(:) !<NODE_LINES(nlp). The local line number of the nlp'th line that is surrounding the node.
    INTEGER(INTG) :: NUMBER_OF_NODE_FACES !<The number of faces surrounding the node in the domain.
    INTEGER(INTG), ALLOCATABLE :: NODE_FACES(:) !<NODE_FACES(nlp). The local face number of the nlp'th face that is surrounding the node.
    LOGICAL :: BOUNDARY_NODE !<Is .TRUE. if the node is on the boundary of the mesh for the domain, .FALSE. if not.
  END TYPE DOMAIN_NODE_TYPE

  !>Contains the topology information for the nodes of a domain
  TYPE DOMAIN_NODES_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this nodes topology information.
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of nodes (excluding ghost nodes) in this domain topology.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of nodes (including ghost nodes) in this domain topology.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_NODES !<The number of global nodes in this domain topology.
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_DERIVATIVES !<The maximum number of derivatives over the nodes in this domain topology.
    TYPE(DOMAIN_NODE_TYPE), POINTER :: NODES(:) !<NODES(np). The pointer to the array of topology information for the nodes of this domain. NODES(np) contains the topological information for the np'th local node of the domain. \todo Change this to allocatable???
    TYPE(TREE_TYPE), POINTER :: NODES_TREE !<A tree mapping the domain local number to the region nodes user number.
  END TYPE DOMAIN_NODES_TYPE

  !>Contains the topology information for a domain
  TYPE DOMAIN_TOPOLOGY_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this topology information.
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES !<The pointer to the topology information for the nodes of this domain.
    TYPE(DOMAIN_DOFS_TYPE), POINTER :: DOFS !<The pointer to the topology information for the dofs of this domain. /todo think about getting rid of domain dofs
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS !<The pointer to the topology information for the elements of this domain.
    TYPE(DOMAIN_FACES_TYPE), POINTER :: FACES !<The pointer to the topology information for the faces of this domain.
    TYPE(DOMAIN_LINES_TYPE), POINTER :: LINES !<The pointer to the topology information for the lines of this domain.
  END TYPE DOMAIN_TOPOLOGY_TYPE

  !
  !================================================================================================================================
  !
  ! Distributed matrix vector types
  !

  !>Contains the information for an adjacent domain for transfering the ghost data of a distributed vector to/from the
  !>current domain.
  TYPE DISTRIBUTED_VECTOR_TRANSFER_TYPE
    TYPE(DISTRIBUTED_VECTOR_CMISS_TYPE), POINTER :: CMISS_VECTOR !<The pointer to the CMISS distributed vector object for this transfer information.
    INTEGER(INTG) :: DATA_TYPE !<The data type of the distributed vector. This is "inherited" from the distributed vector.
    INTEGER(INTG) :: SEND_BUFFER_SIZE !<The size of the buffer to send distributed vector data from the current domain to the adjacent domain.
    INTEGER(INTG) :: RECEIVE_BUFFER_SIZE !<The size of the buffer to receive distributed vector data from the adjacent domain to the current domain.
    INTEGER(INTG) :: SEND_TAG_NUMBER !<The MPI tag number for the data sending from the current domain to the adjacent domain. It is calculated as an offset from the base tag number of the distribued vector.
    INTEGER(INTG) :: RECEIVE_TAG_NUMBER !<The MPI tag number for the data receiving from the adjacent domain to the current domain. It is calculated as an offset from the base tag number of the distribued vector.
    INTEGER(INTG) :: MPI_SEND_REQUEST !<The MPI request pointer for sending data from the current domain to the adjacent domain.
    INTEGER(INTG) :: MPI_RECEIVE_REQUEST !<The MPI request pointer for sending data from the adjacent domain to the current domain.
    INTEGER(INTG), ALLOCATABLE :: SEND_BUFFER_INTG(:) !<The integer buffer for sending the distributed integer vector data from the current domain to the adjacent domain.
    REAL(DP), ALLOCATABLE :: SEND_BUFFER_DP(:) !<The double precision real buffer for sending the distributed real vector data from the current domain to the adjacent domain.
    REAL(SP), ALLOCATABLE :: SEND_BUFFER_SP(:) !<The single precision real buffer for sending the distributed real vector data from the current domain to the adjacent domain.
    LOGICAL, ALLOCATABLE :: SEND_BUFFER_L(:) !<The logical buffer for sending the distributed logical vector data from the current domain to the adjacent domain.
    INTEGER(INTG), ALLOCATABLE :: RECEIVE_BUFFER_INTG(:) !<The integer buffer for receiving the distributed integer vector data from the adjacent domain to the current domain.
    REAL(DP), ALLOCATABLE :: RECEIVE_BUFFER_DP(:) !<The double precision real buffer for receiving the distributed real vector data from the adjacent domain to the current domain.
    REAL(SP), ALLOCATABLE :: RECEIVE_BUFFER_SP(:) !<The single precision real buffer for receiving the distributed real vector data from the adjacent domain to the current domain.
    LOGICAL, ALLOCATABLE :: RECEIVE_BUFFER_L(:) !<The logical buffer for receiving the distributed logical vector data from the adjacent domain to the current domain.  
  END TYPE DISTRIBUTED_VECTOR_TRANSFER_TYPE

  !>Contains information for a CMISS distributed vector
  TYPE DISTRIBUTED_VECTOR_CMISS_TYPE
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR !<A pointer to the distributed vector
    INTEGER(INTG) :: BASE_TAG_NUMBER !<The base number for the MPI tag numbers that will be used to communicate the distributed vector data amongst the domains. The base tag number can be thought of as the identification number for the distributed vector object.
    INTEGER(INTG) :: N !<The size of the distributed vector
    INTEGER(INTG) :: DATA_SIZE !<The size of the distributed vector that is held locally by the domain.
    INTEGER(INTG), ALLOCATABLE :: DATA_INTG(:) !<DATA_INTG(i). The integer data for an integer distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    REAL(DP), ALLOCATABLE :: DATA_DP(:) !<DATA_DP(i). The real data for a double precision real distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    REAL(SP), ALLOCATABLE :: DATA_SP(:) !<DATA_SP(i). The real data for a single precision real distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    LOGICAL, ALLOCATABLE :: DATA_L(:) !<DATA_L(i). The logical data for a logical distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain.  
    TYPE(DISTRIBUTED_VECTOR_TRANSFER_TYPE), ALLOCATABLE :: TRANSFERS(:) !<TRANSFERS(adjacent_domain_idx). The transfer information for the adjacent_domain_idx'th adjacent domain to this domain. 
  END TYPE DISTRIBUTED_VECTOR_CMISS_TYPE

  !>Contains information for a PETSc distributed vector
  TYPE DISTRIBUTED_VECTOR_PETSC_TYPE
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR !<A pointer to the distributed vector
    INTEGER(INTG) :: N !<The number of local components in the vector
    INTEGER(INTG) :: GLOBAL_N !<The number of global components in the vector
    INTEGER(INTG) :: DATA_SIZE  !<The size of the distributed vector that is held locally by the domain.
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_NUMBERS(:) !<GLOBAL_NUMBERS(i). The PETSc global number corresponding to the i'th local number.
    LOGICAL :: USE_OVERRIDE_VECTOR !<Is .TRUE. if the override vector is used instead of the standard vector
    TYPE(PETSC_VEC_TYPE) :: VECTOR !<The PETSc vector
    TYPE(PETSC_VEC_TYPE) :: OVERRIDE_VECTOR !<The PETSc override vector
  END TYPE DISTRIBUTED_VECTOR_PETSC_TYPE
  
  !>Contains the information for a vector that is distributed across a number of domains.
  TYPE DISTRIBUTED_VECTOR_TYPE
    LOGICAL :: VECTOR_FINISHED !<!<Is .TRUE. if the distributed vector has finished being created, .FALSE. if not.
    INTEGER(INTG) :: LIBRARY_TYPE !<The format of the distributed vector \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG) :: GHOSTING_TYPE !<The ghosting type \see DISTRIBUTED_MATRIX_VECTOR_GhostingTypes,DISTRIBUTED_MATRIX_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The pointer for the domain mapping that identifies how the vector is distributed amongst the domains.
    INTEGER(INTG) :: DATA_TYPE !<The type of data for the distributed vector \see DISTRIBUTED_MATRIX_VECTOR_DataTypes 
    TYPE(DISTRIBUTED_VECTOR_CMISS_TYPE), POINTER :: CMISS !<A pointer to the CMISS distributed vector information
    TYPE(DISTRIBUTED_VECTOR_PETSC_TYPE), POINTER :: PETSC !<A pointer to the PETSc distributed vector information
  END TYPE DISTRIBUTED_VECTOR_TYPE

  !>Contains information for a CMISS distributed matrix
  TYPE DISTRIBUTED_MATRIX_CMISS_TYPE
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: DISTRIBUTED_MATRIX !<A pointer to the distributed matrix
    INTEGER(INTG) :: BASE_TAG_NUMBER !<The base number for the MPI tag numbers that will be used to communicate the distributed matrix data amongst the domains. The base tag number can be thought of as the identification number for the distributed matrix object.
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to store the rows corresponding to this domain.
!   TYPE(LINKEDLIST),POINTER :: LIST(:) !< \todo Comment
  END TYPE DISTRIBUTED_MATRIX_CMISS_TYPE

  !>Contains information for a PETSc distributed matrix
  TYPE DISTRIBUTED_MATRIX_PETSC_TYPE
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: DISTRIBUTED_MATRIX !<A pointer to the distributed matrix
    INTEGER(INTG) :: M !<The number of local rows in the PETSc matrix
    INTEGER(INTG) :: N !<The number of local columns in the PETSc matrix
    INTEGER(INTG) :: GLOBAL_M !<The number of global rows in the PETSc matrix
    INTEGER(INTG) :: GLOBAL_N !<The number of global columns in the PETSc matrix
    INTEGER(INTG) :: STORAGE_TYPE !<The storage type (sparsity) of the PETSc matrix
    INTEGER(INTG) :: NUMBER_NON_ZEROS !<The number of non-zeros in the PETSc matrix
    INTEGER(INTG) :: DATA_SIZE !<The size of the allocated data in the PETSc matrix
    INTEGER(INTG) :: MAXIMUM_COLUMN_INDICES_PER_ROW !<The maximum number of column indicies for the rows.
    INTEGER(INTG), ALLOCATABLE :: DIAGONAL_NUMBER_NON_ZEROS(:) !<DIAGONAL_NUMBER_NON_ZEROS(i). The number of non-zeros in the diagonal part of the the i'th row
    INTEGER(INTG), ALLOCATABLE :: OFFDIAGONAL_NUMBER_NON_ZEROS(:) !<OFFDIAGONAL_NUMBER_NON_ZEROS(i). The number of non-zeros in the off diagonal part of the the i'th row
    INTEGER(INTG), ALLOCATABLE :: ROW_INDICES(:) !<ROW_INDICES(i). The row indices for the matrix.
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column indices for the matrix.
    TYPE(LINKEDLIST),POINTER :: LIST(:) !< \todo Comment
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ROW_NUMBERS(:) !<GLOBAL_ROW_NUMBERS(i). The PETSc global row number corresponding to the i'th local row number.
    REAL(DP), POINTER :: DATA_DP(:) !<DATA_DP(i). The real data for the matrix. \todo Is this used???
    LOGICAL :: USE_OVERRIDE_MATRIX !<Is .TRUE. if the override matrix is to be used instead of the standard matrix
    TYPE(PETSC_MAT_TYPE) :: MATRIX !<The PETSc matrix
    TYPE(PETSC_MAT_TYPE) :: OVERRIDE_MATRIX !<The PETSc override matrix
  END TYPE DISTRIBUTED_MATRIX_PETSC_TYPE
  
  !>Contains the information for a matrix that is distributed across a number of domains.
  TYPE DISTRIBUTED_MATRIX_TYPE
    LOGICAL :: MATRIX_FINISHED !<Is .TRUE. if the distributed matrix has finished being created, .FALSE. if not.
    INTEGER(INTG) :: LIBRARY_TYPE !<The library of the distributed matrix \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG) :: GHOSTING_TYPE !<The ghosting type \see DISTRIBUTED_MATRIX_VECTOR_GhostingTypes,DISTRIBUTED_MATRIX_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAPPING !<The pointer for the domain mapping that identifies how the matrix rows are distributed amongst the domains.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOMAIN_MAPPING !<The pointer for the domain mapping that identifies how the matrix columns are distributed amongst the domains.
    INTEGER(INTG) :: DATA_TYPE !<The type of data for the distributed matrix \see DISTRIBUTED_MATRIX_VECTOR_DataTypes
    TYPE(DISTRIBUTED_MATRIX_CMISS_TYPE), POINTER :: CMISS !<A pointer to the CMISS distributed matrix information
    TYPE(DISTRIBUTED_MATRIX_PETSC_TYPE), POINTER :: PETSC !<A pointer to the PETSc distributed matrix information
  END TYPE DISTRIBUTED_MATRIX_TYPE

  !
  !================================================================================================================================
  !
  ! Matrix vector types
  !

  !>Contains information for a vector
  TYPE VECTOR_TYPE
    INTEGER(INTG) :: ID !<The ID of the vector.
    LOGICAL :: VECTOR_FINISHED !<Is .TRUE. if the vector has finished being created, .FALSE. if not.
    INTEGER(INTG) :: N !<The length of the vector
    INTEGER(INTG) :: DATA_TYPE !<The data type of the vector \see MATRIX_VECTOR_DataTypes 
    INTEGER(INTG) :: SIZE !<The size of the data array of the vector
    INTEGER(INTG), ALLOCATABLE :: DATA_INTG(:) !<DATA_INTG(i). The integer data for an integer vector. The i'th component contains the data for the i'th component vector data on the domain. 
    REAL(SP), ALLOCATABLE :: DATA_SP(:) !<DATA_SP(i). The real data for a single precision real vector. The i'th component contains the data for the i'th component vector data on the domain. 
    REAL(DP), ALLOCATABLE :: DATA_DP(:) !<DATA_DP(i). The real data for a double precision real vector. The i'th component contains the data for the i'th component vector data on the domain. 
    LOGICAL, ALLOCATABLE :: DATA_L(:) !<DATA_L(i). The real for a logical vector. The i'th component contains the data for the i'th component vector data on the domain. 
  END TYPE VECTOR_TYPE

  !>Contains information for a matrix
  TYPE MATRIX_TYPE
    INTEGER(INTG) :: ID !<The ID of the matrix
    LOGICAL :: MATRIX_FINISHED !<Is .TRUE. if the matrix has finished being created, .FALSE. if not.
    INTEGER(INTG) :: M !<The number of rows in the matrix
    INTEGER(INTG) :: N !<The number of columns in the matrix
    INTEGER(INTG) :: MAX_M !<The maximum number of columns in the matrix storage
    INTEGER(INTG) :: MAX_N !<The maximum number of rows in the matrix storage
    INTEGER(INTG) :: DATA_TYPE !<The data type of the matrix  \see MATRIX_VECTOR_DataTypes 
    INTEGER(INTG) :: STORAGE_TYPE !<The storage type of the matrix \see MATRIX_VECTOR_StorageTypes 
    INTEGER(INTG) :: NUMBER_NON_ZEROS !<The number of non-zero elements in the matrix 
    INTEGER(INTG) :: SIZE !<The size of the data arrays
    INTEGER(INTG) :: MAXIMUM_COLUMN_INDICES_PER_ROW !<The maximum number of column indicies for the rows.
    INTEGER(INTG), ALLOCATABLE :: ROW_INDICES(:) !<ROW_INDICES(i). The row indices for the matrix storage scheme. \see MATRIX_VECTOR_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column indices for the matrix storage scheme. \see MATRIX_VECTOR_MatrixStorageStructures
    TYPE(LINKEDLIST),POINTER :: LIST(:) !\todo Comment
    INTEGER(INTG), ALLOCATABLE :: DATA_INTG(:) !<DATA_INTG(i). The integer data for an integer matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    REAL(SP), ALLOCATABLE :: DATA_SP(:) !<DATA_SP(i). The real data for a single precision matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    REAL(DP), ALLOCATABLE :: DATA_DP(:) !<DATA_DP(i). The real data for a double precision matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    LOGICAL, ALLOCATABLE :: DATA_L(:) !<DATA_L(i). The logical data for a logical matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
  END TYPE MATRIX_TYPE

  !
  !================================================================================================================================
  !
  ! Domain mapping types
  !
  
  !>Contains the information on an adjacent domain to a domain in a domain mapping. 
  TYPE DOMAIN_ADJACENT_DOMAIN_TYPE
    INTEGER(INTG) :: DOMAIN_NUMBER !<The number of the domain that is adjacent to a domain in a mapping.
    INTEGER(INTG) :: NUMBER_OF_SEND_GHOSTS !<The number of ghost locals in the current domain that are to be sent to this domain.
    INTEGER(INTG) :: NUMBER_OF_RECEIVE_GHOSTS !<The number of ghost locals in the current domain that are to be received from this domain.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_GHOST_SEND_INDICES(:) !<LOCAL_GHOST_SEND_INDICES(i). The local numbers of the ghosts in the current domain that are to be sent to this domain.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_GHOST_RECEIVE_INDICES(:) !<LOCAL_GHOST_RECEIVE_INDICES(i). The local numbers of the ghosts in the current domain that are to be received from this domain.
  END TYPE DOMAIN_ADJACENT_DOMAIN_TYPE
  
  !>Contains the local information for a global mapping number for a domain mapping.
  TYPE DOMAIN_GLOBAL_MAPPING_TYPE
    INTEGER(INTG) :: NUMBER_OF_DOMAINS !<The number of domains that the global number is mapped to a local number in.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NUMBER(:) !<LOCAL_NUMBER(domain_idx). The mapped local number for the domain_idx'th domain for the global number.
    INTEGER(INTG), ALLOCATABLE :: DOMAIN_NUMBER(:) !<DOMAIN_NUMBER(domain_idx). The domain number for the domain_idx'th domain for which the global number is mapped to a local number
    INTEGER(INTG), ALLOCATABLE :: LOCAL_TYPE(:) !<LOCAL_TYPE(domain_idx). The type of local for the domain_idx'th domain for which the global number is mapped to a local number. The types depend on whether the mapped local number in the domain_idx'th domain is an internal, boundary or ghost local number. \see DOMAIN_MAPPINGS_DomainType
  END TYPE DOMAIN_GLOBAL_MAPPING_TYPE

  !>Contains information on the domain mappings (i.e., local and global numberings).
  TYPE DOMAIN_MAPPING_TYPE
    INTEGER(INTG) :: NUMBER_OF_LOCAL !<The number of local numbers in the domain excluding ghost numbers
    INTEGER(INTG) :: TOTAL_NUMBER_OF_LOCAL !<The total number of local numbers in the domain including ghost numbers.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_DOMAIN_LOCAL(:) !<NUMBER_OF_DOMAIN_LOCAL(domain_no). The number of locals for domain_no'th domain. NOTE: the domain_no goes from 0 to the number of domains-1.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_DOMAIN_GHOST(:) !<NUMBER_OF_DOMAIN_GHOST(domain_no). The total number of ghosts for domain_no'th domain. NOTE: the domain_no goes from 0 to the number of domains-1.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL !<The number of global numbers for this mapping.
    INTEGER(INTG) :: NUMBER_OF_DOMAINS !<The number of domains in this mapping.
    INTEGER(INTG) :: NUMBER_OF_INTERNAL !<The number of internal numbers in this mapping.
    INTEGER(INTG) :: NUMBER_OF_BOUNDARY !<The number of boundary numbers in this mapping.
    INTEGER(INTG) :: NUMBER_OF_GHOST !<The number of ghost numbers in this mapping.
    INTEGER(INTG) :: INTERNAL_START !<The start postition in the DOMAIN_LIST for the list of internal numbers
    INTEGER(INTG) :: INTERNAL_FINISH !<The finish postition in the DOMAIN_LIST for the list of internal numbers
    INTEGER(INTG) :: BOUNDARY_START !<The start postition in the DOMAIN_LIST for the list of boundary numbers
    INTEGER(INTG) :: BOUNDARY_FINISH !<The finish postition in the DOMAIN_LIST for the list of boundary numbers
    INTEGER(INTG) :: GHOST_START !<The start postition in the DOMAIN_LIST for the list of ghost numbers
    INTEGER(INTG) :: GHOST_FINISH !<The finish postition in the DOMAIN_LIST for the list of ghost numbers
    INTEGER(INTG), ALLOCATABLE :: DOMAIN_LIST(:) !<DOMAIN_LIST(i). The list of local numbers grouped so that the internal numbers are from INTERNAL_START to INTERNAL_FINISH, the boundary numbers are from BOUNDARY_START to BOUNDARY_FINISH and the ghost numbers are from GHOST_START to GHOST_FINISH
    INTEGER(INTG), ALLOCATABLE :: LOCAL_TO_GLOBAL_MAP(:) !<LOCAL_TO_GLOBAL_MAP(i). The global number for the i'th local number for the mapping.
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), ALLOCATABLE :: GLOBAL_TO_LOCAL_MAP(:) !<GLOBAL_TO_LOCAL_MAP(i). The local information for the i'th global number for the mapping.
    INTEGER(INTG) :: NUMBER_OF_ADJACENT_DOMAINS !<The number of domains that are adjacent to this domain in the mapping.
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_DOMAINS_PTR(:) !<ADJACENT_DOMAINS_PTR(domain_no). The pointer to the list of adjacent domains for domain_no. ADJACENT_DOMAINS_PTR(domain_no) gives the starting position in ADJACENT_DOMAINS_LIST for the first adjacent domain number for domain number domain_no. ADJACENT_DOMAINS_PTR(domain_no+1) gives the last+1 position in ADJACENT_DOMAINS_LIST for the last adjacent domain number for domain number domain_no. NOTE: the index for ADJACENT_DOMAINS_PTR varies from 0 to the number of domains.
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_DOMAINS_LIST(:) !<ADJACENT_DOMAINS_LIST(i). The list of adjacent domains for each domain. The start and end positions for the list for domain number domain_no are given by ADJACENT_DOMAIN_PTR(domain_no) and ADJACENT_DOMAIN_PTR(domain_no+1)-1 respectively.
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), ALLOCATABLE :: ADJACENT_DOMAINS(:) !<ADJACENT_DOMAINS(adjacent_domain_idx). The adjacent domain information for the adjacent_domain_idx'th adjacent domain to this domain. 
  END TYPE DOMAIN_MAPPING_TYPE

  !>Contains information on the domain decomposition mappings.
  TYPE DOMAIN_MAPPINGS_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain decomposition.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS !<Pointer to the element mappings for the domain decomposition.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES !<Pointer to the node mappings for the domain decomposition.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS !<Pointer to the dof mappings for the domain decomposition.
  END TYPE DOMAIN_MAPPINGS_TYPE
  
  !>A pointer to the domain decomposition for this domain.
  TYPE DOMAIN_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the domain decomposition for this domain.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh for this domain.
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The mesh component number of the mesh which this domain was decomposed from.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region that this domain is in. This is "inherited" from the mesh region. 
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS !<The number of dimensions for this domain. This is "inherited" from the mesh.
    INTEGER(INTG), ALLOCATABLE :: NODE_DOMAIN(:) !<NODE_DOMAIN(np). The domain number that the np'th global node is in for the domain decomposition. Note: the domain numbers start at 0 and go up to the NUMBER_OF_DOMAINS-1.
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: MAPPINGS !<Pointer to the mappings for the domain  e.g., global to local and local to global maps
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<Pointer to the topology for the domain.
  END TYPE DOMAIN_TYPE

  !>A buffer type to allow for an array of pointers to a DOMAIN_TYPE.
  TYPE DOMAIN_PTR_TYPE 
    TYPE(DOMAIN_TYPE), POINTER :: PTR !<The pointer to the domain.
  END TYPE DOMAIN_PTR_TYPE

  !
  !================================================================================================================================
  !
  ! Decomposition types
  !
  
  !>Contains the information for a line in a decomposition.
  TYPE DECOMPOSITION_LINE_TYPE
    INTEGER(INTG) :: NUMBER !<The line number in the decomposition.
    INTEGER(INTG) :: XI_DIRECTION !<The Xi direction of the line. Old CMISS name NPL(1,0,nl)
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements that surround (use) this line.
    INTEGER(INTG), ALLOCATABLE :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nel). The local element number of the nel'th element that surrounds (uses) this line. 
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LINES(:) !<ELEMENT_LINES(nel). The local arc number of the nel'th element that surrounds (uses) this line.
    INTEGER(INTG) :: ADJACENT_LINES(0:1) !<ADJACENT_LINES(0:1). The line number of adjacent lines. ADJACENT_LINES(0) is the line number adjacent in the -xi direction. ADJACENT_LINES(1) is the line number adjacent in the +xi direction. Old CMISS name NPL(2..3,0,nl).
    LOGICAL :: BOUNDARY_LINE !<Is .TRUE. if the line is on the boundary of the mesh for the domain, .FALSE. if not.
  END TYPE DECOMPOSITION_LINE_TYPE

  !>Contains the topology information for the lines of a decomposition.
  TYPE DECOMPOSITION_LINES_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this lines topology information.
    INTEGER(INTG) :: NUMBER_OF_LINES !<The number of lines in this decomposition topology.
    TYPE(DECOMPOSITION_LINE_TYPE), ALLOCATABLE :: LINES(:) !<LINES(nl). The pointer to the array of topology information for the lines of this decomposition. LINES(nl) contains the topological information for the nl'th local line of the decomposition.
  END TYPE DECOMPOSITION_LINES_TYPE

  !>Contains the information for a face in a decomposition.
  TYPE DECOMPOSITION_FACE_TYPE
    INTEGER(INTG) :: NUMBER !<The face number in the decomposition.
    INTEGER(INTG) :: XI_DIRECTION !<The Xi direction of the face, the direction of the normal to the face
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements that surround (use) this face.
    INTEGER(INTG), ALLOCATABLE :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nel). The local element number of the nel'th element that surrounds (uses) this face. 
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_FACES(:) !<ELEMENT_FACES(nel). The local arc number of the nel'th element that surrounds (uses) this face.
!    INTEGER(INTG) :: ADJACENT_FACES(0:1) !<ADJACENT_FACES(0:1). The face number of adjacent faces. ADJACENT_FACES(0) is the face number adjacent in the -xi direction. ADJACENT_FACES(1) is the face number adjacent in the +xi direction. Old CMISS name NPL(2..3,0,nl).
    LOGICAL :: BOUNDARY_FACE !<Is .TRUE. if the face is on the boundary of the mesh for the domain, .FALSE. if not.
    INTEGER(INTG) :: ELEMENT_NUMBER !<The element number of the element on which the face is on
  END TYPE DECOMPOSITION_FACE_TYPE

  !>Contains the topology information for the faces of a decomposition.
  TYPE DECOMPOSITION_FACES_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this faces topology information.
    INTEGER(INTG) :: NUMBER_OF_FACES !<The number of faces in this decomposition topology.
    TYPE(DECOMPOSITION_FACE_TYPE), ALLOCATABLE :: FACES(:) !<FACES(nl). The pointer to the array of topology information for the faces of this decomposition. FACES(nl) contains the topological information for the nl'th local face of the decomposition.
  END TYPE DECOMPOSITION_FACES_TYPE

  !>Contains information on the decomposition adjacent elements for a xi coordinate
  TYPE DECOMPOSITION_ADJACENT_ELEMENT_TYPE
    INTEGER(INTG) :: NUMBER_OF_ADJACENT_ELEMENTS !<The number of adjacent elements for the xi coordinate
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:) !<The local element numbers of the elements adjacent to this element for the xi coordinate
  END TYPE DECOMPOSITION_ADJACENT_ELEMENT_TYPE
  
  !>Contains the information for an element in a decomposition.
  TYPE DECOMPOSITION_ELEMENT_TYPE
    INTEGER(INTG) :: LOCAL_NUMBER !<The local element number in the decomposition.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global element number in the mesh of the local element number in the decomposition.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the element.
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS(:) !<ADJACENT_ELEMENTS(-nic:nic). The adjacent elements information in the nic'th xi coordinate direction. Note that -nic gives the adjacent element before the element in the nic'th direction and +nic gives the adjacent element after the element in the nic'th direction. The ni=0 index will give the information on the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    !\todo INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ADJACENT_ELEMENTS(:) !<NUMBER_OF_ADJACENT_ELEMENTS(-ni:ni). The number of elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent element before the element in the ni'th direction and +ni gives the adjacent element after the element in the ni'th direction. The ni=0 index should be 1 for the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    !\todo INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:,:) !<ADJACENT_ELEMENTS(nei,-ni:ni). The local element numbers of the elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent elements before the element in the ni'th direction and +ni gives the adjacent elements after the element in the ni'th direction. The ni=0 index should give the current element number. Old CMISS name NXI(-ni:ni,0:nei,ne).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LINES(:) !<ELEMENT_LINES(nae). The local decomposition line number corresponding to the nae'th local line of the element. Old CMISS name NLL(nae,ne). 
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_FACES(:) !<ELEMENT_FACES(nae). The local decomposition face number corresponding to the nae'th local face of the element. Old CMISS name NLL(nae,ne). 
    LOGICAL :: BOUNDARY_ELEMENT !<Is .TRUE. if the element is on the boundary of the mesh for the domain, .FALSE. if not.
  END TYPE DECOMPOSITION_ELEMENT_TYPE

  !>Contains the topology information for the elements of a decomposition.
  TYPE DECOMPOSITION_ELEMENTS_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this elements topology information.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !<The number of elements excluding ghost elements in this decomposition topology.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ELEMENTS !<The total number of elements in this decomposition topology.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ELEMENTS !<The number of global elements in this decomposition topology.
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: ELEMENTS(:) !<ELEMENTS(ne). The pointer to the array of topology information for the elements of this decomposition. ELEMENTS(ne) contains the topological information for the ne'th local element of the decomposition. \todo Change this to allocatable???
    TYPE(TREE_TYPE), POINTER :: ELEMENTS_TREE !<A tree mapping the decomposition local element number to the decomposition user element number.
  END TYPE DECOMPOSITION_ELEMENTS_TYPE

   !>Contains the topology information for a decomposition
  TYPE DECOMPOSITION_TOPOLOGY_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this topology information.
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS !<The pointer to the topology information for the elements of this decomposition.
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: LINES !<The pointer to the topology information for the lines of this decomposition.
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: FACES !<The pointer to the topology information for the faces of this decomposition.
  END TYPE DECOMPOSITION_TOPOLOGY_TYPE

  !>Contains information on the mesh decomposition. \see OPENCMISS::CMISSDecompositionType
  TYPE DECOMPOSITION_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the domain decomposition. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the domain decomposition in the list of domain decompositions for a particular mesh.
    LOGICAL :: DECOMPOSITION_FINISHED !<Is .TRUE. if the decomposition has finished being created, .FALSE. if not.
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS !<A pointer to the decompositions for this decomposition.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh for this decomposition.
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The component number (index) of the mesh component that this decomposition belongs to (i.e., was generated from).
    INTEGER(INTG) :: DECOMPOSITION_TYPE !<The type of the domain decomposition \see MESH_ROUTINES_DecompositionTypes.
    INTEGER(INTG) :: NUMBER_OF_DOMAINS !<The number of domains that this decomposition contains.
    INTEGER(INTG) :: NUMBER_OF_EDGES_CUT !<For automatically calculated decompositions, the number of edges of the mesh dual graph that were cut for the composition. It provides an indication of the optimally of the automatic decomposition.
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_DOMAIN(:) !<ELEMENT_DOMAIN(ne). The domain number that the ne'th global element is in for the decomposition. Note: the domain numbers start at 0 and go up to the NUMBER_OF_DOMAINS-1.
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the topology for this decomposition.
    TYPE(DOMAIN_PTR_TYPE), POINTER :: DOMAIN(:) !<DOMAIN(mesh_component_idx). A pointer to the domain for mesh component for the domain associated with the computational node. \todo Change this to allocatable???
    LOGICAL :: CALCULATE_FACES !<Boolean flag to determine whether faces should be calculated
    LOGICAL :: CALCULATE_LINES !<Boolean flag to determine whether lines should be calculated
  END TYPE DECOMPOSITION_TYPE

  !>A buffer type to allow for an array of pointers to a DECOMPOSITION_TYPE.
  TYPE DECOMPOSITION_PTR_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: PTR !<The pointer to the domain decomposition. 
  END TYPE DECOMPOSITION_PTR_TYPE

  !>Contains information on the domain decompositions defined on a mesh.
  TYPE DECOMPOSITIONS_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh.
    INTEGER(INTG) :: NUMBER_OF_DECOMPOSITIONS !<The number of decompositions defined on the mesh.
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: DECOMPOSITIONS(:) !<DECOMPOSITIONS(decomposition_idx). The array of pointers to the domain decompositions.
  END TYPE DECOMPOSITIONS_TYPE

  !
  !================================================================================================================================
  !
  ! Field types
  !

  !>Contains information on a physical point in a field.
  TYPE FIELD_PHYSICAL_POINT_TYPE
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FIELD_INTERPOLATED_POINT !<A pointer to the interpolated point of the field that is to be interpolated.
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT !<A pointer to the interpolated point of the geometric field that is to be interpolated.
    INTEGER(INTG) :: PHYSICAL_DERIVATIVE_TYPE !<The type of the physical derivatives that have been interpolated. \see CONSTANTS_PhysicalDerivativeConstants
    REAL(DP), ALLOCATABLE :: VALUES(:) !<VALUES(component_idx). The physical field component values.
  END TYPE FIELD_PHYSICAL_POINT_TYPE

  !>Buffer type to allow for arrays of pointers to FIELD_PHYSICAL_POINT_TYPE
  TYPE FIELD_PHYSICAL_POINT_PTR_TYPE
    TYPE(FIELD_PHYSICAL_POINT_TYPE), POINTER :: PTR
  END TYPE FIELD_PHYSICAL_POINT_PTR_TYPE
  
  !> Contains the interpolated point coordinate metrics. Old CMISS name GL,GU,RG.
  TYPE FIELD_INTERPOLATED_POINT_METRICS_TYPE
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<A pointer to the interpolated point.
    INTEGER(INTG) :: NUMBER_OF_X_DIMENSIONS !<The number of X dimensions.
    INTEGER(INTG) :: NUMBER_OF_XI_DIMENSIONS !<The number of Xi dimensions
    REAL(DP), ALLOCATABLE :: GL(:,:) !<GL(m_xi_idx,n_xi_idx). Covariant metric tensor. Old CMISS name GL.
    REAL(DP), ALLOCATABLE :: GU(:,:) !<GU(m_xi_idx,n_xi_idx). Contravariant metric tensor. Old CMISS name GU.
    REAL(DP), ALLOCATABLE :: DX_DXI(:,:) !<DX_DXI(coord_idx,xi_idx). Rate of change of the X coordinate system wrt the x coordinate system.
    REAL(DP), ALLOCATABLE :: DXI_DX(:,:) !<DXI_DX(xi_idx,coord_idx). Rate of change of the Xi coordinate system wrt the x coordinate system. 
    REAL(DP) :: JACOBIAN !<The Jacobian of the Xi to X coordinate system transformation. Old CMISS name RG.
    INTEGER(INTG) :: JACOBIAN_TYPE !<The type of Jacobian. \see COORDINATE_ROUTINES_JacobianType
  END TYPE FIELD_INTERPOLATED_POINT_METRICS_TYPE

  TYPE FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: PTR
  END TYPE FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE
  
  !>Contains the interpolated value (and the derivatives wrt xi) of a field at a point. Old CMISS name XG.
  TYPE FIELD_INTERPOLATED_POINT_TYPE
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<A pointer to the interpolation parameters of the field that is to be interpolated.
    INTEGER(INTG) :: MAX_PARTIAL_DERIVATIVE_INDEX !<The maximum number of partial derivatives that have been allocated for the values component.
    INTEGER(INTG) :: PARTIAL_DERIVATIVE_TYPE !<The type of the partial derivatives that have been interpolated. PARTIAL_DERIVATIVE_TYPE can be either NO_PART_DERIV, FIRST_PART_DERIV or SECOND_PART_DERIV depending on whether just the field value, the field value and all first derivatives (including cross derivatives) or the first value and all first and second derivatives have been interpolated.
    REAL(DP), ALLOCATABLE :: VALUES(:,:) !<VALUES(component_idx,nu). The interpolated field components and their partial derivatives.
  END TYPE FIELD_INTERPOLATED_POINT_TYPE

  TYPE FIELD_INTERPOLATED_POINT_PTR_TYPE
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: PTR
  END TYPE FIELD_INTERPOLATED_POINT_PTR_TYPE
  
  !>Contains the parameters required to interpolate a field variable within an element. Old CMISS name XE
  TYPE FIELD_INTERPOLATION_PARAMETERS_TYPE
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to be interpolated.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field VARIABLE to be interpolated.
    INTEGER(INTG) :: NUMBER_OF_XI !<The number of xi directions for the interpolation parameters.
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: BASES(:) !<BASES(component_idx). An array to hold a pointer to the basis (if any) used for interpolating the component_idx'th component of the field variable.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_PARAMETERS(:) !<NUMBER_OF_PARAMETERS(component_idx). The number of interpolation parameters used for interpolating the component_idx'th component of the field variable.
    REAL(DP), ALLOCATABLE :: PARAMETERS(:,:) !<PARAMETERS(ns,component_idx). The ns'th interpolation parameter used for interpolating the component_idx'th component of the field variable.
    REAL(DP), ALLOCATABLE :: SCALE_FACTORS(:,:) !<SCALE_FACTORS(ns,component_idx). The scale factors used for scaling the component_idx'th component of the field variable.
  END TYPE FIELD_INTERPOLATION_PARAMETERS_TYPE

  TYPE FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: PTR
  END TYPE FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE
  
  !>Contains the geometric parameters (lines, faces, volumes etc.) for a geometric field decomposition.
  TYPE FIELD_GEOMETRIC_PARAMETERS_TYPE
    INTEGER(INTG) :: NUMBER_OF_LINES !<The number of lines in the field. 
    INTEGER(INTG) :: NUMBER_OF_AREAS !<The number of areas in the field. 
    INTEGER(INTG) :: NUMBER_OF_VOLUMES !<The number of volumes in the field. Inherited from the field decomposition.
    REAL(DP), ALLOCATABLE :: LENGTHS(:) !<LENGTHS(nl). The length of the nl'th line in the field decomposition.
    REAL(DP), ALLOCATABLE :: AREAS(:) !<AREAS(nf). The area of the nf'th face in the field decomposition.
    REAL(DP), ALLOCATABLE :: VOLUMES(:) !<VOLUMES(ne). The volume of the ne'th element in the field decomposition.
    INTEGER(INTG) :: NUMBER_OF_FIELDS_USING !<The number of fields that use these geometric parameters for their scaling. 
    TYPE(FIELD_PTR_TYPE), POINTER :: FIELDS_USING(:) !< FIELDS_USINGS(field_idx). A pointer to the field_idx'th field that uses these geometric parameters for its scaling.
  END TYPE FIELD_GEOMETRIC_PARAMETERS_TYPE

  !>A type to hold the scale factors for the appropriate mesh component of a field. 
  TYPE FIELD_SCALING_TYPE
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The mesh component number of a field variable component that the scaling factors are associated with.
    INTEGER(INTG) :: MAX_NUMBER_OF_DERIVATIVES !<The maximum number of derivatives in the mesh component. 
    INTEGER(INTG) :: MAX_NUMBER_OF_ELEMENT_PARAMETERS !<The maximum number of element parameters in the mesh component. 
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SCALE_FACTORS !<SCALE_FACTORS(nk,np). The scale factor that is applied to the nk'th derivative of the np'th node of the mesh component. \todo  Make scale factors nodally based for now. Will have to revert to element based and extended to be a matrix to allow for a global derivative to be mapped onto many different element derivatives at the points that closes meshes or for inconsistent xi directions
  END TYPE FIELD_SCALING_TYPE

  !>A type to hold the field scalings for the field.
  TYPE FIELD_SCALINGS_TYPE
    INTEGER(INTG) :: SCALING_TYPE !<The type of scaling that is applied to the field. \see FIELD_ROUTINES_ScalingTypes
    INTEGER(INTG) :: NUMBER_OF_SCALING_INDICES !<The number of scaling indices (or sets of scale factors) for the field. In general there will be a set of scale factors (or a scaling index) for each different mesh component that is used by the field variable components.
    TYPE(FIELD_SCALING_TYPE), ALLOCATABLE :: SCALINGS(:) !<SCALINGS(scaling_idx). The scaling factors for the scaling_idx'th set of scaling factors. 
  END TYPE FIELD_SCALINGS_TYPE

  !> A type to hold the mapping from field dof numbers to field parameters (nodes, elements, etc)
  TYPE FIELD_DOF_TO_PARAM_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_DOFS !<The number of degrees-of-freedom for the field.
    INTEGER(INTG), ALLOCATABLE :: DOF_TYPE(:,:) !<DOF_TYPE(i=1..2,ny). The parameter type of the ny'th dof. When i=1 the DOF_TYPE is the type of the dof, i.e., 1=constant field dof, 2=element based field dof, 3=node based field dof, 4=point based field dof. When i=2 the DOF_TYPE gives the nyy'th number of the different field types and is used to index the XXX_DOF2PARAM_MAP arrays.
    INTEGER(INTG) :: NUMBER_OF_CONSTANT_DOFS !<The number of constant degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_ELEMENT_DOFS !<The number of element based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_NODE_DOFS !<The number of node based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_GRID_POINT_DOFS !<The number of grid point based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINT_DOFS !<The number of Gauss point based degrees-of-freedom in the field dofs.
    INTEGER(INTG), ALLOCATABLE :: CONSTANT_DOF2PARAM_MAP(:) !<CONSTANT_DOF2PARAM_MAP(nyy). The mapping from constant field dofs to field parameters for the nyy'th constant field dof. The DOF2PARAM_MAP gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_DOF2PARAM_MAP(:,:) !<ELEMENT_DOF2PARAM_MAP(i=1..2,nyy). The mapping from element based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2PARAM_MAP gives the element number (ne) of the field parameter. When i=2 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.
    INTEGER(INTG), ALLOCATABLE :: NODE_DOF2PARAM_MAP(:,:) !<NODE_DOF2PARAM_MAP(i=1..4,nyy). The mapping from node based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2PARAM_MAP gives the version number (version_idx) of the field parameter. When i=2 the DOF2PARAM_MAP gives the derivative number (derivative_idx) of the field parameter. When i=3 the DOF2PARAM_MAP gives the node number (node_idx) of the field parameter. When i=4 the DOF2PARAM_MAP gives the component number (component_idx) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.
    INTEGER(INTG), ALLOCATABLE :: GRID_POINT_DOF2PARAM_MAP(:,:) !<GRID_POINT_DOF2PARAM_MAP(i=1..2,nyy). The mapping from grid point based field dofs to field parameters for the nyy'th grid point field dof. When i=1 the DOF2PARAM_MAP gives the grid point number (nq) of the field parameter. When i=2 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.  
    INTEGER(INTG), ALLOCATABLE :: GAUSS_POINT_DOF2PARAM_MAP(:,:) !<GAUSS_POINT_DOF2PARAM_MAP(i=1..3,nyy). The mapping from Gauss point based field dofs to field parameters for the nyy'th grid point field dof. When i=1 the DOF2PARAM_MAP gives the Gauss point number (ng) of the field parameter. When i=2 the DOF2PARAM_MAP gives the element number (ne) of the field parameter. When i=3 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.  
  END TYPE FIELD_DOF_TO_PARAM_MAP_TYPE

  !>A type to hold the mapping from a field node derivative's versions to field dof numbers for a particular field variable component.
  TYPE FIELD_NODE_PARAM_TO_DOF_MAP_DERIVATIVE_TYPE
    INTEGER(INTG) :: NUMBER_OF_VERSIONS !<The number of versions for the node derivative parameters of this field variable component.
    INTEGER(INTG), ALLOCATABLE :: VERSIONS(:) !<NODE_PARAM2DOF_MAP%NODES%(node_idx)%DERIVATIVES(derivative_idx)%VERSIONS(version_idx). The field variable dof number of the node_idx'th node's derivative_idx'th derivative's version_idx'th version.
  END TYPE FIELD_NODE_PARAM_TO_DOF_MAP_DERIVATIVE_TYPE

  !>A type to hold the mapping from a field node's derivative to field dof numbers for a particular field variable component.
  TYPE FIELD_NODE_PARAM_TO_DOF_MAP_NODE_TYPE
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES !<The number of derivatives for the node parameter of this field variable component.
    TYPE(FIELD_NODE_PARAM_TO_DOF_MAP_DERIVATIVE_TYPE), ALLOCATABLE :: DERIVATIVES(:) ! The mapping from field node derivative parameter to a dof
  END TYPE FIELD_NODE_PARAM_TO_DOF_MAP_NODE_TYPE

  !>A type to hold the mapping from field nodes to field dof numbers for a particular field variable component.
  TYPE FIELD_NODE_PARAM_TO_DOF_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_NODE_PARAMETERS !<The number of node based field parameters for this field variable component.
    TYPE(FIELD_NODE_PARAM_TO_DOF_MAP_NODE_TYPE), ALLOCATABLE :: NODES(:)  ! The mapping from field node parameter to a dof
  END TYPE FIELD_NODE_PARAM_TO_DOF_MAP_TYPE

  !>A type to hold the mapping from field elements to field dof numbers for a particular field variable component.
  TYPE FIELD_ELEMENT_PARAM_TO_DOF_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_ELEMENT_PARAMETERS !<The number of element based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: ELEMENTS(:) !<ELEMENT_PARAM2DOF_MAP%ELEMENTS(element_idx). The field variable dof number of the element_idx'th element based parameter for this field variable component. \todo Allow for multiple element parameters per element.
  END TYPE FIELD_ELEMENT_PARAM_TO_DOF_MAP_TYPE

  !>A type to hold the mapping from field grid points to field dof numbers for a particular field variable component.
  TYPE FIELD_GRID_POINT_PARAM_TO_DOF_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_GRID_POINT_PARAMETERS !<The number of grid point based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: GRID_POINTS(:) !<GRID_POINT_PARAM2DOF_MAP%GRID_POINTS(grid_point_idx). The field variable dof number of the grid_point_idx'th point based parameter for this field variable component. 
  END TYPE FIELD_GRID_POINT_PARAM_TO_DOF_MAP_TYPE

  !>A type to hold the mapping from field Gauss points to field dof numbers for a particular field variable component.
  TYPE FIELD_GAUSS_POINT_PARAM_TO_DOF_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINT_PARAMETERS !<The number of Gauss point based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: GAUSS_POINTS(:,:) !<GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_point_idx,element_idx). The field variable dof number of the gauss_point_idx'th Gauss point in the element_idx'th element based parameter for this field variable component. 
  END TYPE FIELD_GAUSS_POINT_PARAM_TO_DOF_MAP_TYPE

  !>A type to hold the mapping from field parameters (nodes, elements, etc) to field dof numbers for a particular field variable component.
  TYPE FIELD_PARAM_TO_DOF_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_CONSTANT_PARAMETERS !<The number of constant field parameters for this field variable component. Note: this is currently always 1 but is included for completeness and to allow for multiple constants per field variable component in the future.
    INTEGER(INTG) :: CONSTANT_PARAM2DOF_MAP !<The field variable dof number of the constant parameter for this field variable component.
    TYPE(FIELD_ELEMENT_PARAM_TO_DOF_MAP_TYPE) :: ELEMENT_PARAM2DOF_MAP !> A type to hold the mapping from field element parameters to field dof numbers
    TYPE(FIELD_NODE_PARAM_TO_DOF_MAP_TYPE) :: NODE_PARAM2DOF_MAP !> A type to hold the mapping from field node parameters to field dof numbers
    TYPE(FIELD_GRID_POINT_PARAM_TO_DOF_MAP_TYPE) :: GRID_POINT_PARAM2DOF_MAP !> A type to hold the mapping from grid point element parameters to field dof numbers
    TYPE(FIELD_GAUSS_POINT_PARAM_TO_DOF_MAP_TYPE) :: GAUSS_POINT_PARAM2DOF_MAP !> A type to hold the mapping from field gauss point parameters to field dof numbers
  END TYPE FIELD_PARAM_TO_DOF_MAP_TYPE

  !>Contains information for a component of a field variable.
  TYPE FIELD_VARIABLE_COMPONENT_TYPE
    INTEGER(INTG) :: COMPONENT_NUMBER !<The number of the field variable component.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable for this component.
    TYPE(VARYING_STRING) :: COMPONENT_LABEL !<The label for the field variable component
    INTEGER(INTG) :: INTERPOLATION_TYPE !<The interpolation type of the field variable component \see FIELD_ROUTINES_InterpolationTypes
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The mesh component of the field decomposition for this field variable component.
    INTEGER(INTG) :: SCALING_INDEX !<The index into the defined field scalings for this field variable component.
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain of the field decomposition for this field variable component.
    INTEGER(INTG) :: MAX_NUMBER_OF_INTERPOLATION_PARAMETERS !<The maximum number of interpolations parameters in an element for a field variable component.
    TYPE(FIELD_PARAM_TO_DOF_MAP_TYPE) :: PARAM_TO_DOF_MAP !<The mapping of the field parameters to the field dofs for this field variable component.
  END TYPE FIELD_VARIABLE_COMPONENT_TYPE

  !>A type to hold the parameter sets for a field.
  TYPE FIELD_PARAMETER_SET_TYPE
    INTEGER(INTG) :: SET_INDEX !<The global set index (from 1 to the TYPES::FIELD_PARAMETER_SETS_TYPE::NUMBER_OF_PARAMETER_SETS) that this parameter set corresponds to.
    INTEGER(INTG) :: SET_TYPE !<The user set type (index) (from 1 to FIELD_ROUTINES::FIELD_NUMBER_OF_SET_TYPES) that this parameter set \see FIELD_ROUTINES_ParameterSetTypes
  !###      corresponds to.
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: PARAMETERS !<A pointer to the distributed vector that contains the field parameters for this field parameter set.
  END TYPE FIELD_PARAMETER_SET_TYPE
  
  !>A buffer type to allow for an array of pointers to a FIELD_PARAMETER_SET_TYPE.
  TYPE FIELD_PARAMETER_SET_PTR_TYPE
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PTR !<The pointer to the field parameter set. 
  END TYPE FIELD_PARAMETER_SET_PTR_TYPE

  !>A type to store the parameter sets for a field.
  TYPE FIELD_PARAMETER_SETS_TYPE    
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable that these parameter sets are defined on.
    INTEGER(INTG) :: NUMBER_OF_PARAMETER_SETS !<The number of parameter sets that are currently defined on the field.
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE), POINTER :: SET_TYPE(:) !<SET_TYPE(set_type_idx). A pointer to an array of pointers to the field set types. SET_TYPE(set_type_idx)%PTR is a pointer to the parameter set type for the set_type_idx'th parameter set. set_type_idx can vary from 1 to FIELD_ROUTINES::FIELD_NUMBER_OF_SET_TYPES. The value of the pointer will be NULL if the parameter set corresponding to the set_type_idx'th parameter set has not yet been created for the field.
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE), POINTER :: PARAMETER_SETS(:) !<PARAMETER_SETS(set_type_idx). A pointer to an array of pointers to the parameter sets that have been created on the field. PARAMETER_SET(set_type_idx)%PTR is a pointer to the parameter set type for the set_type_idx'th parameter set that has been created. set_type_idx can vary from 1 to the number of parameter set types that have currently been created for the field i.e., TYPES::FIELD_PARAMETER_SETS_TYPE::NUMBER_OF_PARAMETER_SETS.
  END TYPE FIELD_PARAMETER_SETS_TYPE
  
  !>Contains information for a field variable defined on a field.
  TYPE FIELD_VARIABLE_TYPE
    INTEGER(INTG) :: VARIABLE_NUMBER !<The number of the field variable
    INTEGER(INTG) :: VARIABLE_TYPE !<The type of the field variable. \see FIELD_ROUTINES_VariableTypes
    TYPE(VARYING_STRING) :: VARIABLE_LABEL !<The label for the variable
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field for this field variable.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region for this field variable.
    INTEGER(INTG) :: DIMENSION !<The dimension of the field variable. \see FIELD_ROUTINES_DimensionTypes
    INTEGER(INTG) :: DATA_TYPE !<The data type of the field variable.  \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
    INTEGER(INTG) :: DOF_ORDER_TYPE !<The order of the DOF's in the field variable \see FIELD_ROUTINES_DOFOrderTypes,FIELD_ROUTINES
    INTEGER(INTG) :: MAX_NUMBER_OF_INTERPOLATION_PARAMETERS !<The maximum number of interpolation parameters in an element for a field variable. 
    INTEGER(INTG) :: NUMBER_OF_DOFS !<Number of local degress of freedom for this field variable (excluding ghosted dofs). Old CMISS name NYNR(0,0,nc,nr,nx).
    INTEGER(INTG) :: TOTAL_NUMBER_OF_DOFS !<Number of local degrees of freedom for this field variable (including ghosted dofs). Old CMISS name NYNR(0,0,nc,nr,nx).
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_DOFS !<Number of global degrees of freedom for this field variable. Old CMISS name NYNR(0,0,nc,nr,nx).
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<Domain mapping for this variable. 
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !<The number of components in the field variable.
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), ALLOCATABLE :: COMPONENTS(:) !<COMPONENTS(component_idx). The array of field variable components.
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE) :: DOF_TO_PARAM_MAP !<The mappings for the field dofs to the field parameters
    TYPE(FIELD_PARAMETER_SETS_TYPE) :: PARAMETER_SETS !<The parameter sets for the field variable
  END TYPE FIELD_VARIABLE_TYPE
  
  !>A buffer type to allow for an array of pointers to a FIELD_VARIABLE_TYPE.
  TYPE FIELD_VARIABLE_PTR_TYPE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: PTR !<The pointer to the field variable. 
  END TYPE FIELD_VARIABLE_PTR_TYPE

  !>A type to temporarily hold (cache) the user modifiable values which are used to create a field. 
  TYPE FIELD_CREATE_VALUES_CACHE_TYPE
    LOGICAL :: LABEL_LOCKED !<Is .TRUE. if the field label has been locked, .FALSE. if not.
    LOGICAL :: DECOMPOSITION_LOCKED !<Is .TRUE. if the field decomposition has been locked, .FALSE. if not.
    LOGICAL :: DEPENDENT_TYPE_LOCKED !<Is .TRUE. if the field dependent type has been locked, .FALSE. if not.
    LOGICAL :: NUMBER_OF_VARIABLES_LOCKED !<Is .TRUE. if the number of field variables has been locked, .FALSE. if not.
    LOGICAL :: GEOMETRIC_FIELD_LOCKED !<Is .TRUE. if the geometric field has been locked, .FALSE. if not.
    LOGICAL :: SCALING_TYPE_LOCKED !<Is .TRUE. if the scaling type has been locked, .FALSE. if not.        
    LOGICAL :: TYPE_LOCKED !<Is .TRUE. if the field type has been locked, .FALSE. if not.        
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). The cache of the variable type for the given variable_idx of the field. \see FIELD_ROUTINES_VariableTypes
    LOGICAL :: VARIABLE_TYPES_LOCKED !<Is .TRUE. if the variable types have been locked, .FALSE. if not.
    TYPE(VARYING_STRING), ALLOCATABLE :: VARIABLE_LABELS(:) !<VARIABLE_LABELS(variable_type_idx). The variable label for the variable_type_idx'th variable type of the field.
    LOGICAL, ALLOCATABLE :: VARIABLE_LABELS_LOCKED(:) !<VARIABLE_LABELS_LOCKED(variable_type_idx). Is .TRUE. if the variable label for the variable_type_idx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: DIMENSION(:) !<DIMENSION(variable_type_idx). The cache of the variable dimension for the variable_type_idx'th variable type of the field. \see FIELD_ROUTINES_DimensionTypes
    LOGICAL, ALLOCATABLE :: DIMENSION_LOCKED(:) !<DIMENSION_LOCKED(variable_type_idx). Is .TRUE. if the dimension for the variable_type_idx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: DATA_TYPES(:) !<DATA_TYPES(variable_type_idx). The cache of the variable data type for the variable_type_idx'th variable type of the field. \see FIELD_ROUTINES_DataTypes
    LOGICAL, ALLOCATABLE :: DATA_TYPES_LOCKED(:) !<DATA_TYPES_LOCKED(variable_type_idx). Is .TRUE. if the data type for the variable_type_idx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: DOF_ORDER_TYPES(:) !<DOF_ORDER_TYPES(variable_type_idx). The cache of the variable dof order type for the variable_type_idx'th variable type of the field. \see FIELD_ROUTINES_DataTypes
    LOGICAL, ALLOCATABLE :: DOF_ORDER_TYPES_LOCKED(:) !<DOF_ORDER_TYPES_LOCKED(variable_type_idx). Is .TRUE. if the dof order type for the variable_type_idx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_COMPONENTS(:) !<NUMBER_OF_COMPONENTS(variable_type_idx). The number of components in the field for the variable_type_idx'th field variable type.
    LOGICAL, ALLOCATABLE :: NUMBER_OF_COMPONENTS_LOCKED(:) !<NUMBER_OF_COMPONENTS_LOCKED(variable_type_idx). Is .TRUE. if the number of components has been locked for the variable_type_idx'th variable type, .FALSE. if not.
    TYPE(VARYING_STRING), ALLOCATABLE :: COMPONENT_LABELS(:,:) !<COMPONENT_LABELS(component_idx,variable_type_idx). The cache of the component label for the given component and variable type of the field.
    LOGICAL, ALLOCATABLE :: COMPONENT_LABELS_LOCKED(:,:) !<COMPONENT_LABELS_LOCKED(component_idx,variable_type_idx). Is .TRUE. if the component label of the component_idx'th component of the variable_type_idx'th varible type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: INTERPOLATION_TYPE(:,:) !<INTERPOLATION_TYPES(component_idx,variable_type_idx). The cache of the interpolation type for the given component and variable type of the field. \see FIELD_ROUTINES_InterpolationTypes
    LOGICAL, ALLOCATABLE :: INTERPOLATION_TYPE_LOCKED(:,:) !<INTERPOLATION_TYPES(component_idx,variable_type_idx). Is .TRUE. if the interpolation type of the component_idx'th component of the variable_type_idx'th varible type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENT_NUMBER(:,:) !<MESH_COMPONENT_NUMBER(component_idx,varaible_type_idx). The cache of the mesh component number for the given component and variable type of the field.
    LOGICAL, ALLOCATABLE :: MESH_COMPONENT_NUMBER_LOCKED(:,:) !<MESH_COMPONENT_NUMBER_LOCKED(component_idx,variable_type_idx). Is .TRUE. if the mesh component number of the component_idx'th component of the variable_type_idx'th varible type has been locked, .FALSE. if not.
  END TYPE FIELD_CREATE_VALUES_CACHE_TYPE

  !>Contains information for a field defined on a region. \see OPENCMISS::CMISSFieldType
  TYPE FIELD_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the field in the list of fields for a region.
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the field. The user number must be unique.
    TYPE(VARYING_STRING) :: LABEL !<The label for the field
    LOGICAL :: FIELD_FINISHED !<Is .TRUE. if the field has finished being created, .FALSE. if not.
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointer to the fields for this region.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the field. If the field are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE!<A pointer to the interface containing the field. If the field are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: TYPE !<The type of the field. \see FIELD_ROUTINES_FieldTypes
    INTEGER(INTG) :: DEPENDENT_TYPE !<The dependent type of the field. \see FIELD_ROUTINES_DependentTypes
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition of the mesh for which the field is defined on.
    INTEGER(INTG) :: NUMBER_OF_VARIABLES !<The number of variable types in the field. Old CMISS name NCT(nr,nx)
    TYPE(FIELD_VARIABLE_PTR_TYPE), ALLOCATABLE :: VARIABLE_TYPE_MAP(:) !<VARIABLE_TYPE_MAP(variable_idx). The map from the available field variable types to the field variable types that are defined for the field. variable_idx varies from 1 to FIELD_ROUTINES::FIELD_NUMBER_OF_VARIABLE_TYPES. If the particular field variable type has not been defined on the field then the VARIABLE_TYPE_MAP will be NULL. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_VARIABLE_TYPE), ALLOCATABLE :: VARIABLES(:) !<VARIABLES(variable_idx). The array of field variables. 
    TYPE(FIELD_SCALINGS_TYPE) :: SCALINGS !<The scaling parameters for the field
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field that this field uses. If the field itself is a geometric field then this will be a pointer back to itself.
    TYPE(FIELD_GEOMETRIC_PARAMETERS_TYPE), POINTER :: GEOMETRIC_FIELD_PARAMETERS !<If the field is a geometric field the pointer to the geometric parameters (lines, areas, volumes etc.). If the field is not a geometric field the pointer is NULL.
    TYPE(FIELD_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<The create values cache for the field.
  END TYPE FIELD_TYPE

  !>A buffer type to allow for an array of pointers to a FIELD_TYPE.
  TYPE FIELD_PTR_TYPE
    TYPE(FIELD_TYPE), POINTER :: PTR !<The pointer to the field.  
  END TYPE FIELD_PTR_TYPE

  !>Contains information on the fields defined on a region.
  TYPE FIELDS_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the fields. If the fields are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the fields. If the fields are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: NUMBER_OF_FIELDS !<The number of fields defined on the region.
    TYPE(FIELD_PTR_TYPE), POINTER :: FIELDS(:) !<FIELDS(fields_idx). The array of pointers to the fields.
  END TYPE FIELDS_TYPE

  !
  !================================================================================================================================
  !
  ! Equations matrices types
  !
  
  !>Contains information for an element matrix.
  TYPE ELEMENT_MATRIX_TYPE
    INTEGER(INTG) :: EQUATIONS_MATRIX_NUMBER !<The equations matrix number that this element matrix belongs to.
    INTEGER(INTG) :: STRUCTURE_TYPE !<The structure type of the element matrix. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The current number of rows in the element matrix.
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The current number of columns in the element matrix.
    INTEGER(INTG) :: MAX_NUMBER_OF_ROWS !<The maximum (allocated) number of rows in the element matrix.
    INTEGER(INTG) :: MAX_NUMBER_OF_COLUMNS !<The maximu (allocated) number of columns in the element matrix.
    INTEGER(INTG), ALLOCATABLE :: ROW_DOFS(:) !<ROW_DOFS(i). The equations row that the i'th row of the element matrix belongs to.
    INTEGER(INTG), ALLOCATABLE :: COLUMN_DOFS(:) !<COLUMN_DOFS(j). The equations column that the j'th column of the element matrix bleongs to.
    REAL(DP), ALLOCATABLE :: MATRIX(:,:) !<MATRIX(i,j). The vlaue of the i'th row and the j'th column of the element matrix.
  END TYPE ELEMENT_MATRIX_TYPE

  !>Contains information for an element vector.
  TYPE ELEMENT_VECTOR_TYPE
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The current number of rows in the element vector
    INTEGER(INTG) :: MAX_NUMBER_OF_ROWS !<The maximum (allocated) number of rows in the element vecotr
    INTEGER(INTG), ALLOCATABLE :: ROW_DOFS(:) !<ROW_DOFS(i). The equations row that the i'th row of the element vector belongs to
    REAL(DP), ALLOCATABLE :: VECTOR(:) !<VECTOR(i). The value of the i'th row of the element vector
  END TYPE ELEMENT_VECTOR_TYPE

  !>Contains information for an nodal matrix.
  TYPE NodalMatrixType
    INTEGER(INTG) :: equationsMatrixNumber !<The equations matrix number that this nodal matrix belongs to.
    INTEGER(INTG) :: structureType !<The structure type of the nodal matrix. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
    INTEGER(INTG) :: numberOfRows !<The current number of rows in the nodal matrix.
    INTEGER(INTG) :: numberOfColumns !<The current number of columns in the nodal matrix.
    INTEGER(INTG) :: maxNumberOfRows !<The maximum (allocated) number of rows in the nodal matrix.
    INTEGER(INTG) :: maxNumberOfColumns !<The maximum (allocated) number of columns in the nodal matrix.
    INTEGER(INTG), ALLOCATABLE :: rowDofs(:) !<rowDofx(i). The equations row that the i'th row of the nodal matrix belongs to.
    INTEGER(INTG), ALLOCATABLE :: columnDofs(:) !<columnDofs(j). The equations column that the j'th column of the nodal matrix bleongs to.
    REAL(DP), ALLOCATABLE :: matrix(:,:) !<matrix(i,j). The vlaue of the i'th row and the j'th column of the nodal matrix.
  END TYPE NodalMatrixType

  !>Contains information for an nodal vector.
  TYPE NodalVectorType
    INTEGER(INTG) :: numberOfRows !<The current number of rows in the nodal vector
    INTEGER(INTG) :: maxNumberOfRows !<The maximum (allocated) number of rows in the nodal vecotr
    INTEGER(INTG), ALLOCATABLE :: rowDofs(:) !<rowDofs(i). The equations row that the i'th row of the nodal vector belongs to
    REAL(DP), ALLOCATABLE :: vector(:) !<vector(i). The value of the i'th row of the nodal vector
  END TYPE NodalVectorType

  !>Contains information about an equations matrix.
  TYPE EQUATIONS_MATRIX_TYPE
    INTEGER(INTG) :: MATRIX_NUMBER !<The number of the equations matrix
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES !<A pointer to the dynamic equations matrices for the dynamic equation matrix.
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES !<A pointer to the linear equations matrices for the linear equation matrix.
    INTEGER(INTG) :: STORAGE_TYPE !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: STRUCTURE_TYPE !<The structure (sparsity) type for this matrix
    LOGICAL :: LUMPED !<Is .TRUE. if the equations matrix is lumped
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in this equations matrix
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if this equations matrix is to be updated
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this equations matrix has not been assembled
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the distributed equations matrix data
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix for this equations matrix
    TYPE(NodalMatrixType) :: NodalMatrix !<The nodal matrix for this equations matrix
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: TEMP_VECTOR !<Temporary vector used for assembly. 
  END TYPE EQUATIONS_MATRIX_TYPE

  !>A buffer type to allow for an array of pointers to a EQUATIONS_MATRIX_TYPE \see TYPES::EQUATIONS_MATRIX_TYPE.
  TYPE EQUATIONS_MATRIX_PTR_TYPE
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: PTR !<A pointer to the equations matrix.
  END TYPE EQUATIONS_MATRIX_PTR_TYPE
 
  !>Contains information on the Jacobian matrix for nonlinear problems
  TYPE EQUATIONS_JACOBIAN_TYPE
    INTEGER(INTG) :: JACOBIAN_NUMBER !<The equations Jacobian matrix number
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer back to the nonlinear matrices for this Jacobian
    INTEGER(INTG) :: STORAGE_TYPE !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: STRUCTURE_TYPE !<The structure (sparsity) type for this matrix
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in this global matrix
    LOGICAL :: UPDATE_JACOBIAN !<Is .TRUE. if this Jacobian matrix is to be updated
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: JACOBIAN !<A pointer to the distributed jacobian matrix data
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this Jacobian matrix has not been assembled
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_JACOBIAN !<The element matrix for this Jacobian matrix. This is not used if the Jacobian is not supplied.
    TYPE(NodalMatrixType) :: NodalJacobian !<The nodal matrix for this Jacobian matrix. This is not used if the Jacobian is not supplied.
    INTEGER(INTG) :: JACOBIAN_CALCULATION_TYPE !<The calculation type (analytic of finite difference) of the Jacobian.
  END TYPE EQUATIONS_JACOBIAN_TYPE

  !>A buffer type to allow for an array of pointers to a EQUATIONS_JACOBIAN_TYPE \see TYPES::EQUATIONS_JACOBIAN_TYPE.
  TYPE EQUATIONS_JACOBIAN_PTR_TYPE
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: PTR
  END TYPE EQUATIONS_JACOBIAN_PTR_TYPE

  !>Contains information of the dynamic matrices for equations matrices
  TYPE EQUATIONS_MATRICES_DYNAMIC_TYPE
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer back to the equations matrices.
    INTEGER(INTG) :: NUMBER_OF_DYNAMIC_MATRICES !<The number of dynamic equations matrices defined for the equations set.
    TYPE(EQUATIONS_MATRIX_PTR_TYPE), ALLOCATABLE :: MATRICES(:) !<MATRICES(matrix_idx)%PTR contains the information on the matrix_idx'th dynamic equations matrix.
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: TEMP_VECTOR !<Temporary vector used for assembly. 
  END TYPE EQUATIONS_MATRICES_DYNAMIC_TYPE

  !>Contains information of the linear matrices for equations matrices
  TYPE EQUATIONS_MATRICES_LINEAR_TYPE
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer back to the equations matrices.
    INTEGER(INTG) :: NUMBER_OF_LINEAR_MATRICES !<The number of linear equations matrices defined for the equations set.
    TYPE(EQUATIONS_MATRIX_PTR_TYPE), ALLOCATABLE :: MATRICES(:) !<MATRICES(matrix_idx)%PTR contains the information on the matrix_idx'th linear equations matrix.
  END TYPE EQUATIONS_MATRICES_LINEAR_TYPE

  !>Contains information of the nolinear matrices and vectors for equations matrices
  TYPE EQUATIONS_MATRICES_NONLINEAR_TYPE
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer back to the equations matrices.
    INTEGER(INTG) :: NUMBER_OF_JACOBIANS !<The number of Jacobian matrices for the equations set.
    TYPE(EQUATIONS_JACOBIAN_PTR_TYPE), ALLOCATABLE :: JACOBIANS(:) !<JACOBIANS(matrix_idx)%PTR is a pointer to the matrix_idx'th Jacobian matrix for nonlinear equations
    LOGICAL :: UPDATE_RESIDUAL !<Is .TRUE. if the equtions residual vector is to be updated
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this residual vector has not been assembled
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL !<A pointer to the distributed residual vector for nonlinear equations
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_RESIDUAL !<The element residual information for nonlinear equations. Old CMISS name RE1
    TYPE(NodalVectorType) :: NodalResidual !<The nodal residual information for nonlinear equations.
    INTEGER(INTG) :: NodalResidualCalculated !<The number of the nodal the residual is calculated for, or zero if it isn't calculated
    INTEGER(INTG) :: ELEMENT_RESIDUAL_CALCULATED !<The number of the element the residual is calculated for, or zero if it isn't calculated
  END TYPE EQUATIONS_MATRICES_NONLINEAR_TYPE

  !>Contains information of the RHS vector for equations matrices
  TYPE EQUATIONS_MATRICES_RHS_TYPE
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer back to the equations matrices.
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the equtions rhs vector is to be updated
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this rhs vector has not been assembled
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the distributed global rhs vector data \todo rename this RHS_VECTOR
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element rhs information
    TYPE(NodalVectorType) :: NodalVector !<The nodal rhs information
  END TYPE EQUATIONS_MATRICES_RHS_TYPE
  
  !>Contains information of the source vector for equations matrices
  TYPE EQUATIONS_MATRICES_SOURCE_TYPE
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer back to the equations matrices.
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the equtions rhs vector is to be updated
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this source vector has not been assembled
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the distributed source vector data \todo rename this SOURCE_VECTOR
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element source information
    TYPE(NodalVectorType) :: NodalVector !<The nodal source information
  END TYPE EQUATIONS_MATRICES_SOURCE_TYPE
  
  !>Contains information on the equations matrices and vectors
  TYPE EQUATIONS_MATRICES_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer back to the equations
    LOGICAL :: EQUATIONS_MATRICES_FINISHED !<Is .TRUE. if the equations matrices have finished being created, .FALSE. if not.
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping for the equations matrices.
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping for the equations matrices
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of local rows (excluding ghost rows) in the distributed equations matrices and vectors
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ROWS !<The number of local rows (including ghost rows) in the distributed equations matrices and vectors
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ROWS !<The number of global rows in the distributed equations matrices and vectors
    !Equations matrices components
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES !<A pointer to the dynamic matrices information for the equations matrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES !<A pointer to the linear matrices information for the equations matrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer to the nonlinear matrices and vectors information for the equations matrices
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR !<A pointer to the RHS vector information for the equations matrices
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR !<A pointer to the source vector information for the equations matrices
  END TYPE EQUATIONS_MATRICES_TYPE

  !
  !================================================================================================================================
  !
  ! Equations mapping types
  !
  
   !>Contains the information about the mapping of a variable DOF to an equations matrix column
  TYPE VAR_TO_EQUATIONS_COLUMN_MAP_TYPE
    INTEGER(INTG), ALLOCATABLE :: COLUMN_DOF(:) !<COLUMN_DOF(dof_idx). The equations column number for this equations matrix that the dof_idx'th variable DOF is mapped to.  
  END TYPE VAR_TO_EQUATIONS_COLUMN_MAP_TYPE

  !>Contains the mapping for a dependent variable type to the equations matrices
  TYPE VAR_TO_EQUATIONS_MATRICES_MAP_TYPE
    INTEGER(INTG) :: VARIABLE_INDEX !<The variable index for this variable to equations matrices map
    INTEGER(INTG) :: VARIABLE_TYPE !<The variable type for this variable to equations matrices map
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the field variable for this variable to equations matrices map
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS_MATRICES !<The number of equations matrices (linear or dynamic) this variable type is mapped to. If the number is -1 the variable is mapped to the RHS vector. If the number is zero then this variable type is not involved in the equations set and the rest of the type is not allocated.
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_MATRIX_NUMBERS(:) !<EQUATIONS_MATRIX_NUMBERS(i). The equations matrix number for the i'th matrix that this variable type is mapped to.
    TYPE(VAR_TO_EQUATIONS_COLUMN_MAP_TYPE), ALLOCATABLE :: DOF_TO_COLUMNS_MAPS(:) !<DOF_TO_COLUMNS_MAPS(i). The variable dof to equations columns for the i'th equations matrix.
    INTEGER(INTG), ALLOCATABLE :: DOF_TO_ROWS_MAP(:) !<DOF_TO_ROWS_MAP(dof_idx). The row number that the dof_idx'th variable dof is mapped to.
  END TYPE VAR_TO_EQUATIONS_MATRICES_MAP_TYPE

  !>Contains information for mapping an equations matrix to a field variable.
  TYPE EQUATIONS_MATRIX_TO_VAR_MAP_TYPE
    INTEGER(INTG) :: MATRIX_NUMBER !<The equations matrix number
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX !<A pointer to the equations matrix
    INTEGER(INTG) :: VARIABLE_TYPE !<The dependent variable type mapped to this equations matrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the field variable that is mapped to this equations matrix
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in this equations matrix.
    REAL(DP) :: MATRIX_COEFFICIENT !<The multiplicative coefficent for the matrix in the equation set
    INTEGER(INTG), ALLOCATABLE :: COLUMN_TO_DOF_MAP(:) !<COLUMN_TO_DOF_MAP(column_idx). The variable DOF that the column_idx'th column of this equations matrix is mapped to.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOFS_MAPPING !<A pointer to the column dofs domain mapping for the matrix variable
  END TYPE EQUATIONS_MATRIX_TO_VAR_MAP_TYPE

  !>Contains information for mapping field variables to the dynamic matrices in the equations set of the mapping
  TYPE EQUATIONS_MAPPING_DYNAMIC_TYPE
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG) :: NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES !<The number of dynamic equations matrices in this mapping

    INTEGER(INTG) :: STIFFNESS_MATRIX_NUMBER !<The matrix number of the dynamic stiffness matrix. 0 if there is no dynamic stiffness matrix
    INTEGER(INTG) :: DAMPING_MATRIX_NUMBER !<The matrix number of the dynamic damping matrix. 0 if there is no dynamic damping matrix
    INTEGER(INTG) :: MASS_MATRIX_NUMBER !<The matrix number of the dynamic mass matrix. 0 if there is no dynamic mass matrix
    INTEGER(INTG) :: DYNAMIC_VARIABLE_TYPE !<The variable type involved in the equations matrix mapping.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DYNAMIC_VARIABLE !<A pointer to the variable that is mapped to the dynamic matrices.
!!TODO: just make this the size of the number of matrix variables used (i.e. 1) rather than the field number of variable types???
    TYPE(VAR_TO_EQUATIONS_MATRICES_MAP_TYPE), ALLOCATABLE :: VAR_TO_EQUATIONS_MATRICES_MAPS(:) !<VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type_idx). The equations matrices mapping for the variable_type_idx'th variable type.
    TYPE(EQUATIONS_MATRIX_TO_VAR_MAP_TYPE), ALLOCATABLE :: EQUATIONS_MATRIX_TO_VAR_MAPS(:) !<EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx). The mappings for the matrix_idx'th equations matrix.
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(:) !<EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(row_idx). The row mappings for the row_idx'th row of the equations matrices to the dynamic variable.
  END TYPE EQUATIONS_MAPPING_DYNAMIC_TYPE

  !>Contains information for mapping field variables to the linear matrices in the equations set of the mapping
  TYPE EQUATIONS_MAPPING_LINEAR_TYPE
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG) :: NUMBER_OF_LINEAR_EQUATIONS_MATRICES !<The number of linear equations matrices in this mapping
    INTEGER(INTG) :: NUMBER_OF_LINEAR_MATRIX_VARIABLES !<The number of dependent variables involved in the linear equations matrix mapping
    INTEGER(INTG), ALLOCATABLE :: LINEAR_MATRIX_VARIABLE_TYPES(:) !<LINEAR_MATRIX_VARIABLE_TYPES(i). The variable type of the i'th variable type involved in the equations linear matrix mapping.
!!TODO: just make this the size of the number of matrix variables rather than the field number of variable types and merge matrix variable types above???
    TYPE(VAR_TO_EQUATIONS_MATRICES_MAP_TYPE), ALLOCATABLE :: VAR_TO_EQUATIONS_MATRICES_MAPS(:) !<VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type_idx). The equations matrices mapping for the variable_type_idx'th variable type.
    TYPE(EQUATIONS_MATRIX_TO_VAR_MAP_TYPE), ALLOCATABLE :: EQUATIONS_MATRIX_TO_VAR_MAPS(:) !<EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx). The mappings for the matrix_idx'th equations matrix.
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(:,:) !<EQUATIONS_ROW_TO_VARIABLE_DOFS_MAPS(row_idx,variable_type_idx). The row mappings for the row_idx'th row of the equations matrices to the variable_type_idx'th variable.
  END TYPE EQUATIONS_MAPPING_LINEAR_TYPE

  !>Contains the mapping from the Jacobian back to the nonlinear residual variables.
  TYPE EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE
    INTEGER(INTG) :: JACOBIAN_NUMBER !<The equations Jacobian matrix number
    INTEGER(INTG) :: VARIABLE_TYPE !<The dependent variable type mapped to this equations matrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the field variable that is mapped to this equations matrix
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN !<A pointer to the equations matrix for this variable
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in this equations matrix.
    REAL(DP) :: JACOBIAN_COEFFICIENT !<The multiplicative coefficent for the matrix in the equation set
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP(:) !<COLUMN_TO_DOF_MAP(column_idx). The variable DOF that the column_idx'th column of this equations matrix is mapped to.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOFS_MAPPING !<A pointer to the column dofs domain mapping for the matrix variable
  END TYPE EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE

  !>Contains the mapping for a dependent variable type to the nonlinear Jacobian matrix
  TYPE VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE
    INTEGER(INTG) :: JACOBIAN_NUMBER !<The equations Jacobian matrix number
    INTEGER(INTG) :: VARIABLE_TYPE !<The variable type for this variable to equations matrices map
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the field variable for this variable to equations matrices map
    INTEGER(INTG), ALLOCATABLE :: DOF_TO_COLUMNS_MAP(:) !<DOF_TO_COLUMNS_MAP(dof_idx). The Jacobian column number for dof_idx'th variable dof
    INTEGER(INTG), ALLOCATABLE :: DOF_TO_ROWS_MAP(:) !<DOF_TO_ROWS_MAP(dof_idx). The row number that the dof_idx'th variable dof is mapped to.
  END TYPE VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE

  !>Contains information on the equations mapping for nonlinear matrices i.e., how a field variable is mapped to residual
  !>vectors, and how the field variables are mapped to the rows and columns of the associated Jacobian matrices of the equations set
  !>of this equations mapping.
  !>There may be multiple residual variables with a Jacobian matrix for each variable
  TYPE EQUATIONS_MAPPING_NONLINEAR_TYPE
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG) :: NUMBER_OF_RESIDUAL_VARIABLES !<The number of residual variables in this mapping
    TYPE(FIELD_VARIABLE_PTR_TYPE), ALLOCATABLE :: RESIDUAL_VARIABLES(:) !<RESIDUAL_VARIABLES(variable_idx). The variable_idx'th residual variable.
    TYPE(VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE), ALLOCATABLE :: VAR_TO_JACOBIAN_MAP(:) !<VAR_TO_JACOBIAN_MAP(variable_idx). The mapping from the residual variable to the Jacobain matrix for the variable_idx'th residual variable.
    TYPE(EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE), ALLOCATABLE :: JACOBIAN_TO_VAR_MAP(:) !<JACOBIAN_TO_VAR_MAP(jacobian_idx). The mapping from the Jacobian matrix to the residual variables for the jacobian_idx'th Jacobian.
    REAL(DP) :: RESIDUAL_COEFFICIENT !<The multiplicative coefficient applied to the residual vector
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP(:) !<EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP(row_idx). The mapping from the row_idx'th row of the equations to the source dof.
  END TYPE EQUATIONS_MAPPING_NONLINEAR_TYPE

  !>Contains information on the equations mapping for a RHS i.e., how a field variable is mapped to the RHS vector for
  !>the equations set of this equations mapping.
  TYPE EQUATIONS_MAPPING_RHS_TYPE
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG) :: RHS_VARIABLE_TYPE !<The variable type number mapped to the RHS vector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: RHS_VARIABLE !<A pointer to the variable that is mapped to the RHS vector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: RHS_VARIABLE_MAPPING !<A pointer to the RHS variable domain mapping
    REAL(DP) :: RHS_COEFFICIENT !<The multiplicative coefficient applied to the RHS vector
    INTEGER(INTG), ALLOCATABLE :: RHS_DOF_TO_EQUATIONS_ROW_MAP(:) !<RHS_DOF_TO_EQUATIONS_ROW_MAP(residual_dof_idx). The mapping from the rhs_dof_idx'th RHS dof in the rhs variable to the equations row.   
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_ROW_TO_RHS_DOF_MAP(:) !<EQUATIONS_ROW_TO_RHS_DOF_MAP(row_idx). The mapping from the row_idx'th row of the equations to the RHS dof.   
  END TYPE EQUATIONS_MAPPING_RHS_TYPE

  !>Contains information on the equations mapping for a source i.e., how a field variable is mapped to the source vector for
  !>the equations set of this equations mapping.
  TYPE EQUATIONS_MAPPING_SOURCE_TYPE
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG) :: SOURCE_VARIABLE_TYPE !<The variable type number mapped from the source vector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: SOURCE_VARIABLE !<A pointer to the source variable 
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: SOURCE_VARIABLE_MAPPING !<A pointer to the domain mapping for the source variable.
    REAL(DP) :: SOURCE_COEFFICIENT !<The multiplicative coefficient applied to the source vector
    INTEGER(INTG), ALLOCATABLE :: SOURCE_DOF_TO_EQUATIONS_ROW_MAP(:) !<SOURCE_DOF_TO_EQUATIONS_ROW_MAP(source_dof_idx). The mapping from the source_dof_idx'th source dof in the source variable to the equations row.   
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_ROW_TO_SOURCE_DOF_MAP(:) !<EQUATIONS_ROW_TO_SOURCE_DOF_MAP(row_idx). The mapping from the row_idx'th row of the equations to the source dof.
  END TYPE EQUATIONS_MAPPING_SOURCE_TYPE

  !>Contains information on the create values cache for the equations mapping. Because we do not want to allocate and deallocate
  !>large data structures as the equations mapping options are changed between create start and create finish we cache the
  !>important information and the allocate and process the data structures at create finish.
  TYPE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE
    INTEGER(INTG) :: NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES !<The number of dynamic matrices in the equations mapping
    INTEGER(INTG) :: DYNAMIC_STIFFNESS_MATRIX_NUMBER !<The dynamic matrix number corresponding to the dynamic stiffness matrix
    INTEGER(INTG) :: DYNAMIC_DAMPING_MATRIX_NUMBER !<The dynamic matrix number corresponding to the dynamic damping matrix
    INTEGER(INTG) :: DYNAMIC_MASS_MATRIX_NUMBER !<The dynamic matrix number corresponding to the dynamic mass matrix
    INTEGER(INTG) :: DYNAMIC_VARIABLE_TYPE !<The dependent variable type mapped to the dynamic equations matrices.
    REAL(DP), ALLOCATABLE :: DYNAMIC_MATRIX_COEFFICIENTS(:) !<DYNAMIC_MATRIX_COEFFICIENTS(matrix_idx). The coefficient of the matrix_idx'th dynamic matrix in the equations set.
    INTEGER(INTG) :: NUMBER_OF_LINEAR_EQUATIONS_MATRICES !<The number of linear matrices in the equations mapping
    INTEGER(INTG), ALLOCATABLE :: LINEAR_MATRIX_VARIABLE_TYPES(:) !<LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx). The dependent variable type mapped to the matrix_idx'th linear equations matrix.
    REAL(DP), ALLOCATABLE :: LINEAR_MATRIX_COEFFICIENTS(:) !<LINEAR_MATRIX_COEFFICIENTS(matrix_idx). The coefficient of the matrix_idx'th linear matrix in the equations set.
    INTEGER(INTG) :: NUMBER_OF_RESIDUAL_VARIABLES !<The number of residual variables for the nonlinear equations.
    INTEGER(INTG), ALLOCATABLE :: RESIDUAL_VARIABLE_TYPES(:) !<RESIDUAL_VARIABLE_TYPES(jacobian_idx). The type of the jacobian_idx'th Jacobian.
    REAL(DP) :: RESIDUAL_COEFFICIENT !<The coefficient multiplying the residual vector.
    INTEGER(INTG) :: RHS_VARIABLE_TYPE !<The dependent variable type mapped to the rhs vector.
    REAL(DP) :: RHS_COEFFICIENT !<The coefficient multiplying the RHS vector.
    INTEGER(INTG) :: SOURCE_VARIABLE_TYPE !<The source variable type mapped to the source vector
    REAL(DP) :: SOURCE_COEFFICIENT !<The coefficient multiplying the source vector.
  END TYPE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE

  !>Contains information on the equations mapping i.e., how field variable DOFS are mapped to the rows and columns of a number
  !>of equations matrices.
  TYPE EQUATIONS_MAPPING_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations for this equations mapping
    LOGICAL :: EQUATIONS_MAPPING_FINISHED !<Is .TRUE. if the equations mapping has finished being created, .FALSE. if not.
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices associated with this equations mapping.
    !Row mappings
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of local rows (excluding ghost rows) in the equations matrices
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ROWS !<The number of local rows (including ghost rows) in the equations matrices
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ROWS !<The number of global rows in the equations matrices
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOFS_MAPPING !<The domain mapping for the equations rows
    !Equations mapping components
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING !<A pointer to the equations mapping for dynamic matrices
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING !<A pointer to the equations mapping for the linear matrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING !<A pointer to the equations mapping for the nonlinear matrices and vectors
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING !<A pointer to the equations mapping for the RHS vector
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING !<A pointer to the equations mapping for the source vector
    !Create values cache
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<The create values cache for the equations mapping
  END TYPE EQUATIONS_MAPPING_TYPE

  !
  !================================================================================================================================
  !
  ! Equations types
  !
  
  !>Contains information on the interpolation for the equations
  TYPE EQUATIONS_INTERPOLATION_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field for the equations.
    TYPE(FIELD_TYPE), POINTER :: FIBRE_FIELD !<A pointer to the fibre field for the equations (if one is defined).
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<A pointer to the dependent field for the equations 
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD !<A pointer to the independent field for the equations 
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD !<A pointer to the material field for the equations (if one is defined).
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD !<A pointer to the source field for the equations (if one is defined).
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: GEOMETRIC_INTERP_PARAMETERS(:) !<GEOMETRIC_INTERP_PARAMETERS(field_variable_type). A pointer to the field_variable_type'th geometric interpolation parameters for the equations.
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: FIBRE_INTERP_PARAMETERS(:) !<FIBRE_INTERP_PARAMETERS(field_variable_type). A pointer to the fibre interpolation parameters for the equations (if a fibre field is defined). 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: DEPENDENT_INTERP_PARAMETERS(:) !<DEPENDENT_INTERP_PARAMETERS(field_variable_type). A pointer to the field_variable_type'th dependent interpolation parameters for the equations. 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INDEPENDENT_INTERP_PARAMETERS(:) !<INDEPENDENT_INTERP_PARAMETERS(field_variable_type). A pointer to the field_variable_type'th independent interpolation parameters for the equations. 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: MATERIALS_INTERP_PARAMETERS(:) !<MATERIALS_INTERP_PARAMETERS(field_variable_type). A pointer to the field_variable_type'th material interpolation parameters for the equations (if a material field is defined). 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: SOURCE_INTERP_PARAMETERS(:) !<SOURCE_INTERP_PARAMETERS(field_variable_type). A pointer to the field_variable_type'th source interpolation parameters for the equations (if a source field is defined). 
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: GEOMETRIC_INTERP_POINT(:) !<GEOMETRIC_INTERP_POINT(field_variable_type). A pointer to the field_variable_type'th geometric interpolated point information for the equations. 
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: FIBRE_INTERP_POINT(:) !<FIBRE_INTERP_POINT(field_variable_type). A pointer to the field_variable_type'th fibre interpolated point information for the equations (if a fibre field is defined). 
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: DEPENDENT_INTERP_POINT(:) !<DEPENDENT_INTERP_POINT(field_variable_type). A pointer to the field_variable_type'th dependent interpolated point information for the equations. 
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INDEPENDENT_INTERP_POINT(:) !<INDEPENDENT_INTERP_POINT(field_variable_type). A pointer to the field_variable_type'th independent interpolated point information for the equations. 
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: MATERIALS_INTERP_POINT(:) !<MATERIALS_INTERP_POINT(field_variable_type). A pointer to the field_variable_type'th material interpolated point information for the equations (if a material field is defined). 
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: SOURCE_INTERP_POINT(:) !<SOURCE_INTERP_POINT(field_variable_type). A pointer to the field_variable_type'th source interpolated point information for the equations (if a source field is defined).
    TYPE(FIELD_PHYSICAL_POINT_PTR_TYPE), POINTER :: DEPENDENT_PHYSICAL_POINT(:) !<DEPENDENT_PHYSICAL_POINT(field_variable_type). A pointer to the field_variable_type'th dependent physical point information for the equations. 
     TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: DEPENDENT_INTERP_POINT_METRICS(:) !<DEPENDENT_INTERP_POINT_METRICS(field_variable_type). A pointer to the field_variable_type'th dependent interpolated point metrics information 
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: INDEPENDENT_INTERP_POINT_METRICS(:) !<INDEPENDENT_INTERP_POINT_METRICS(field_variable_type). A pointer to the field_variable_type'th independent interpolated point metrics information 
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: GEOMETRIC_INTERP_POINT_METRICS(:) !<GEOMETRIC_INTERP_POINT_METRICS(field_variable_type). A pointer to the field_variable_type'th geometric interpolated point metrics information 
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: FIBRE_INTERP_POINT_METRICS(:) !<FIBRE_INTERP_POINT_METRICS(field_variable_type). A pointer to the field_variable_type'th fibre interpolated point metrics information 
  END TYPE EQUATIONS_INTERPOLATION_TYPE

  !>Contains information about the equations in an equations set. \see OPENCMISS::CMISSEquationsType
  TYPE EQUATIONS_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations_set
    LOGICAL :: EQUATIONS_FINISHED !<Is .TRUE. if the equations have finished being created, .FALSE. if not.
    INTEGER(INTG) :: LINEARITY !<The equations linearity type \see EQUATIONS_ROUTINES_LinearityTypes,EQUATIONS_ROUTINES
    INTEGER(INTG) :: TIME_DEPENDENCE !<The equations time dependence type \see EQUATIONS_ROUTINES_TimeDepedenceTypes,EQUATIONS_ROUTINES
    INTEGER(INTG) :: OUTPUT_TYPE !<The output type for the equations \see EQUATIONS_ROUTINES_EquationsOutputTypes,EQUATIONS_ROUTINES
    INTEGER(INTG) :: SPARSITY_TYPE !<The sparsity type for the equation matrices of the equations \see EQUATIONS_ROUTINES_EquationsSparsityTypes,EQUATIONS_ROUTINES
    INTEGER(INTG) :: LUMPING_TYPE !<The lumping type for the equation matrices of the equations \see EQUATIONS_ROUTINES_EquationsLumpingTypes,EQUATIONS_ROUTINES
    TYPE(EQUATIONS_INTERPOLATION_TYPE), POINTER :: INTERPOLATION !<A pointer to the interpolation information used in the equations.
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping for the equations.
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices and vectors used for the equations.
  END TYPE EQUATIONS_TYPE

  TYPE EQUATIONS_PTR_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: PTR
  END TYPE EQUATIONS_PTR_TYPE

  !
  !================================================================================================================================
  !
  ! Boundary conditions types
  !

  !>Contains information on the boundary conditions for a dependent field variable
  TYPE BOUNDARY_CONDITIONS_VARIABLE_TYPE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions for this boundary conditions variable
    INTEGER(INTG) :: VARIABLE_TYPE !<The type of variable for this variable boundary conditions
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the field variable for this boundary condition variable
    INTEGER(INTG), ALLOCATABLE :: DOF_TYPES(:) !<DOF_TYPES(dof_idx). The general boundary condition type (eg. fixed or free) of the dof_idx'th dof in the dependent field variable. \see OPENCMISS_BoundaryConditionsTypes,OPENCMISS
    INTEGER(INTG), ALLOCATABLE :: CONDITION_TYPES(:) !<CONDITION_TYPES(dof_idx). The specific boundary condition type (eg. incremented pressure) of the dof_idx'th dof of the dependent field variable, which might be specific to an equation set. The solver routines should not need to use this array, and should only need the DOF_TYPES array. \see OPENCMISS_BoundaryConditionsDOFTypes,OPENCMISS
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: DIRICHLET_BOUNDARY_CONDITIONS  !<A pointer to the dirichlet boundary condition type for this boundary condition variable
    INTEGER(INTG) :: NUMBER_OF_DIRICHLET_CONDITIONS !<Stores the number of dirichlet conditions associated with this variable
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions
    TYPE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE), POINTER :: PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS !<A pointer to the pressure incremented condition type for this boundary condition variable
    INTEGER(INTG), ALLOCATABLE :: DOF_COUNTS(:) !<DOF_COUNTS(CONDITION_TYPE): The number of DOFs that have a CONDITION_TYPE boundary condition set
    LOGICAL, ALLOCATABLE :: parameterSetRequired(:) !<parameterSetRequired(PARAMETER_SET) is true if any boundary condition has been set that requires the PARAMETER_SET field parameter set
  END TYPE BOUNDARY_CONDITIONS_VARIABLE_TYPE

  !>A buffer type to allow for an array of pointers to a VARIABLE_BOUNDARY_CONDITIONS_TYPE \see TYPES::VARIABLE_BOUNDARY_CONDITIONS_TYPE
  TYPE BOUNDARY_CONDITIONS_VARIABLE_PTR_TYPE
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: PTR !<A pointer to the boundary conditions variable
  END TYPE BOUNDARY_CONDITIONS_VARIABLE_PTR_TYPE

  !>Contains information on the boundary conditions for the solver equations. \see OPENCMISS::CMISSBoundaryConditionsType
  TYPE BOUNDARY_CONDITIONS_TYPE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations.
    LOGICAL :: BOUNDARY_CONDITIONS_FINISHED !<Is .TRUE. if the boundary conditions for the equations set has finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES !<The number of boundary conditions variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_PTR_TYPE), ALLOCATABLE :: BOUNDARY_CONDITIONS_VARIABLES(:) !<BOUNDARY_CONDITIONS_VARIABLES(variable_idx). BOUNDARY_CONDITIONS_VARIABLES(variable_idx)%PTR is the pointer to the variable_idx'th boundary conditions variable. variable_idx ranges from 1 to NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
    INTEGER(INTG) :: neumannMatrixSparsity !<The sparsity type of the Neumann integration matrices. \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  END TYPE BOUNDARY_CONDITIONS_TYPE

  !>A buffer type to allow for an array of pointers to a BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE \see TYPES::BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE
  TYPE BOUNDARY_CONDITIONS_SPARSITY_INDICES_PTR_TYPE
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: PTR !<A pointer to the boundary conditions sparsity indices type
  END TYPE BOUNDARY_CONDITIONS_SPARSITY_INDICES_PTR_TYPE

  !> Contains information on dofs with associated dirichlet conditions and corresponding non-zero elements in the equations matrices
  TYPE BOUNDARY_CONDITIONS_DIRICHLET_TYPE
    INTEGER(INTG), ALLOCATABLE :: DIRICHLET_DOF_INDICES(:)  !<DIRICHLET_DOF_INDICES(idx). Stores the dof_idx of the dofs which are subject to a dirichlet boundary condition \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_PTR_TYPE), ALLOCATABLE :: LINEAR_SPARSITY_INDICES(:,:) !<LINEAR_SPARSITY_INDICES(equ_set_idx,equ_matrix_idx). Stores the indices of the non-zero elements of the equ_set_idx'th equation set and equ_matrix_idx'th linear equation matrix in the columns corresponding to the dofs which are subject to a dirichlet boundary condition
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_PTR_TYPE), ALLOCATABLE :: DYNAMIC_SPARSITY_INDICES(:,:) !<DYNAMIC_SPARSITY_INDICES(equ_set_idx,equ_matrix_idx). Stores the indices of the non-zero elements of the equ_set_idx'th equation set and equ_matrix_idx'th dynamic equation matrix in the columns corresponding to the dofs which are subject to a dirichlet boundary condition
  END TYPE BOUNDARY_CONDITIONS_DIRICHLET_TYPE

  !> Contains information on indices of non-zero elements with associated dirichlet conditions
  !> Indices stored in compressed column format without a values array
  TYPE BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE
    INTEGER(INTG), ALLOCATABLE :: SPARSE_ROW_INDICES(:) !<SPARSE_ROW_INDICES(SPARSE_COLUMN_INDICES(column_idx)). Between SPARSE_COLUMN_INDICES(column_idx) and SPARSE_COLUMN_INDICES(column_idx+1)-1 are the row indices of non-zero elements of the 'column_idx'th column
    INTEGER(INTG), ALLOCATABLE :: SPARSE_COLUMN_INDICES(:) !<SPARSE_COLUMN_INDICES(column_idx). Between SPARSE_COLUMN_INDICES(column_idx) and SPARSE_COLUMN_INDICES(column_idx+1)-1 are the row indices of non-zero elements of the 'column_idx'th column
  END TYPE BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE

  !>Contains information used to integrate Neumann boundary conditions
  TYPE BoundaryConditionsNeumannType
    INTEGER(INTG), ALLOCATABLE :: setDofs(:) !<setDofs(neumann_idx): the global dof for the neumann_idx'th Neumann condition
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: integrationMatrix !<The N matrix that multiples the point values vector q to give the integrated values f. Number of rows equals number of local dofs, and number of columns equals number of set point DOFs.
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: pointValues !<The vector of set point values q
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: pointDofMapping !<The domain mapping for DOFs with Neumann point conditions set.
  END TYPE BoundaryConditionsNeumannType

  !>Contains information on dofs associated with pressure incremented conditions
  TYPE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE
    INTEGER(INTG), ALLOCATABLE :: PRESSURE_INCREMENTED_DOF_INDICES(:)  !<PRESSURE_INCREMENTED_DOF_INDICES(idx). Stores the dof_idx of the dofs which are subject to a pressure incremented boundary condition \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  END TYPE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE

  !
  !================================================================================================================================
  !
  ! Equations set types
  !

  !>Contains information on the setup information for an equations set
  TYPE EQUATIONS_SET_SETUP_TYPE
    INTEGER(INTG) :: SETUP_TYPE !<The setup type for the equations set setup \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG) :: ACTION_TYPE !<The action type for the equations set setup \see EQUATIONS_SET_CONSTANTS_SetupActionTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG) :: FIELD_USER_NUMBER !<The user number for the field for the equations set setup.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field for the equations set setup.
    INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE !<The analytic function type to use.
  END TYPE EQUATIONS_SET_SETUP_TYPE  
  
  !>Contains information on the geometry for an equations set
  TYPE EQUATIONS_SET_GEOMETRY_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<The geometric field for this equations set.
    TYPE(FIELD_TYPE), POINTER :: FIBRE_FIELD !<The fibre field for this equations set if one is defined. If no fibre field is defined the pointer is NULL.
  END TYPE EQUATIONS_SET_GEOMETRY_TYPE

  TYPE EQUATIONS_SET_MATERIALS_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    LOGICAL :: MATERIALS_FINISHED !<Is .TRUE. if the materials for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: MATERIALS_FIELD_AUTO_CREATED !<Is .TRUE. if the materials field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD !<A pointer to the materials field for the equations set if one is defined. If no material field is defined the pointer is NULL.
  END TYPE EQUATIONS_SET_MATERIALS_TYPE

  !>Contains information on the dependent variables for the equations set.
  TYPE EQUATIONS_SET_DEPENDENT_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    LOGICAL :: DEPENDENT_FINISHED !<Is .TRUE. if the dependent variables for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: DEPENDENT_FIELD_AUTO_CREATED !<Is .TRUE. if the dependent field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<A pointer to the dependent field for the equations set.
  END TYPE EQUATIONS_SET_DEPENDENT_TYPE

  !>Contains information on the independent variables for the equations set.
  TYPE EQUATIONS_SET_INDEPENDENT_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    LOGICAL :: INDEPENDENT_FINISHED !<Is .TRUE. if the independent variables for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: INDEPENDENT_FIELD_AUTO_CREATED !<Is .TRUE. if the independent field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD !<A pointer to the independent field for the equations set.
  END TYPE EQUATIONS_SET_INDEPENDENT_TYPE

  !>Contains information on the source for the equations set.
  TYPE EQUATIONS_SET_SOURCE_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    LOGICAL :: SOURCE_FINISHED !<Is .TRUE. if the source for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: SOURCE_FIELD_AUTO_CREATED !<Is .TRUE. if the source field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD !<A pointer to the source field for the equations set if one is defined. If no source is defined the pointer is NULL.
  END TYPE EQUATIONS_SET_SOURCE_TYPE
  
  !>Contains information on the analytic setup for the equations set.
  TYPE EQUATIONS_SET_ANALYTIC_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE !<The analytic function identifier
    LOGICAL :: ANALYTIC_FINISHED !<Is .TRUE. if the analytic setup for the problem has finished being created, .FALSE. if not.
    LOGICAL :: ANALYTIC_FIELD_AUTO_CREATED !<Is .TRUE. if the analytic field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD !<A pointer to the analytic field for the equations set if one is defined. If no source is defined the pointer is NULL.
    REAL(DP) :: ANALYTIC_TIME !<The time value to use for analytic evaluations.
    REAL(DP) :: ANALYTIC_USER_PARAMS(20)  !<A small array that can be used to hold various parameters often required in analytic problems. \todo should this be allocated?
  END TYPE EQUATIONS_SET_ANALYTIC_TYPE

  TYPE EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set.
    LOGICAL :: EQUATIONS_SET_FIELD_FINISHED !<Is .TRUE. if the equations set field for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: EQUATIONS_SET_FIELD_AUTO_CREATED !<Is .TRUE. if the equations set field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD !<A pointer to the equations set field for the equations set.
  END TYPE EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE

  !>Contains information on an equations set. \see OPENCMISS::CMISSEquationsSetType
  TYPE EQUATIONS_SET_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user identifying number of the equations set
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global index of the equations set in the region.
    LOGICAL :: EQUATIONS_SET_FINISHED !<Is .TRUE. if the equations set have finished being created, .FALSE. if not.
    TYPE(EQUATIONS_SETS_TYPE), POINTER :: EQUATIONS_SETS !<A pointer back to the equations sets
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer back to the region containing the equations set.
    
    INTEGER(INTG) :: CLASS !<The equations set specification class identifier
    INTEGER(INTG) :: TYPE !<The equations set specification type identifier
    INTEGER(INTG) :: SUBTYPE !<The equations set specification subtype identifier
    
    INTEGER(INTG) :: SOLUTION_METHOD !<The solution method for the equations set \see EQUATIONS_ROUTINES_SolutionMethods 
    
    TYPE(EQUATIONS_SET_GEOMETRY_TYPE) :: GEOMETRY !<The geometry information for the equations set.
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: MATERIALS !<A pointer to the materials information for the equations set.
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: SOURCE !<A pointer to the source information for the equations set.
    TYPE(EQUATIONS_SET_DEPENDENT_TYPE) :: DEPENDENT !<The depedent variable information for the equations set.
    TYPE(EQUATIONS_SET_INDEPENDENT_TYPE), POINTER :: INDEPENDENT !<A pointer to the indepedent field information for the equations set.
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: ANALYTIC !<A pointer to the analytic setup information for the equations set.
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations information for the equations set
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary condition information for the equations set.
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE) :: EQUATIONS_SET_FIELD !<A pointer to the equations set field for the equations set.
  END TYPE EQUATIONS_SET_TYPE
  
  !>A buffer type to allow for an array of pointers to a EQUATIONS_SET_TYPE \see TYPES::EQUATIONS_SET_TYPE
  TYPE EQUATIONS_SET_PTR_TYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: PTR !<A pointer to the equations set.
  END TYPE EQUATIONS_SET_PTR_TYPE
  
  TYPE EQUATIONS_SETS_TYPE
    TYPE(REGION_TYPE) , POINTER :: REGION !<A pointer to the region containing the equations sets
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS_SETS !<The number of equations sets
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: EQUATIONS_SETS(:) !<EQUATIONS_SETS(equations_set_idx). EQUATIONS_SETS(equations_set_idx)%PTR is the pointer to the equations_set_idx'th equations set.
  END TYPE EQUATIONS_SETS_TYPE

  !
  !================================================================================================================================
  !
  ! Interface types
  
  !>Contains information about an interface matrix.
  TYPE INTERFACE_MATRIX_TYPE
    INTEGER(INTG) :: MATRIX_NUMBER !<The number of the interface matrix
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices for the interface matrix.
    INTEGER(INTG) :: STORAGE_TYPE !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: STRUCTURE_TYPE !<The structure (sparsity) type for this matrix
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of rows in this interface matrix
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ROWS !<The number of rows in this interface matrix
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if this interface matrix is to be updated
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this interface matrix has not been assembled
    LOGICAL :: HAS_TRANSPOSE !<Is .TRUE. if this interface matrix has has transpose
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the distributed interface matrix data
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: MATRIX_TRANSPOSE !<A pointer to the distributed interface matrix transpose data
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix for this interface matrix
  END TYPE INTERFACE_MATRIX_TYPE

  !>A buffer type to allow for an array of pointers to a INTERFACE_MATRIX_TYPE \see TYPES::INTERFACE_MATRIX_TYPE.
  TYPE INTERFACE_MATRIX_PTR_TYPE
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: PTR !<A pointer to the interface matrix.
  END TYPE INTERFACE_MATRIX_PTR_TYPE

  !>Contains information of the RHS vector for interface matrices
  TYPE INTERFACE_RHS_TYPE
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer back to the interface matrices.
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the interface rhs vector is to be updated
    LOGICAL :: FIRST_ASSEMBLY !<Is .TRUE. if this rhs vector has not been assembled
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RHS_VECTOR !<A pointer to the distributed global rhs vector data 
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element rhs information
  END TYPE INTERFACE_RHS_TYPE
  
  !>Contains information on the interface matrices 
  TYPE INTERFACE_MATRICES_TYPE
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer back to the interface equations
    LOGICAL :: INTERFACE_MATRICES_FINISHED !<Is .TRUE. if the interface  matrices have finished being created, .FALSE. if not.
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface equations mapping for the interface equations matrices.
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping for the interface equations matrices
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of local columns in the interface matrices
    INTEGER(INTG) :: TOTAL_NUMBER_OF_COLUMNS !<The total number of local columns in the interface matrices
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_COLUMNS !<The number of global columns in the interface matrices
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_MATRICES !<The number of interfaces matrices defined for the interface condition.
    TYPE(INTERFACE_MATRIX_PTR_TYPE), ALLOCATABLE :: MATRICES(:) !<MATRICES(matrix_idx)%PTR contains the information on the matrix_idx'th  interface matrix.
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR !<A pointer to the RHS vector information for the interface matrices.
  END TYPE INTERFACE_MATRICES_TYPE

  !>Contains information on interface variable mapping for an interface matrix.
  TYPE INTERFACE_MATRIX_TO_VAR_MAP_TYPE
    INTEGER(INTG) :: MATRIX_NUMBER !<The interface matrix number
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set containing the dependent variable that is mapped to this interface matrix.
    INTEGER(INTG) :: VARIABLE_TYPE !<The dependent variable type mapped to this interface matrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the field variable that is mapped to this interface matrix
    INTEGER(INTG) :: MESH_INDEX !<The mesh index for the matrix in the interface.
    REAL(DP) :: MATRIX_COEFFICIENT !<The multiplicative coefficent for the matrix.    
    LOGICAL :: HAS_TRANSPOSE !<.TRUE. if the interface matrix has a tranpose, .FALSE. if not.  
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of rows  in this interface matrix.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ROWS !<The total number of rows in this interface matrix.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ROWS !<The global number of rows in this interface matrix.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOFS_MAPPING !<A pointer to the domain mapping for the row dofs.
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_DOF_TO_ROW_MAP(:) !<VARIABLE_DOF_TO_ROW_MAP(dof_idx). The mapping from the dof_idx'th variable dof to the rows on the interface matrix
  END TYPE INTERFACE_MATRIX_TO_VAR_MAP_TYPE

  TYPE INTERFACE_MAPPING_RHS_TYPE
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer back to the interface mapping
    INTEGER(INTG) :: RHS_VARIABLE_TYPE !<The variable type number mapped to the RHS vector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: RHS_VARIABLE !<A pointer to the variable that is mapped to the RHS vector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: RHS_VARIABLE_MAPPING !<A pointer to the RHS variable domain mapping
    REAL(DP) :: RHS_COEFFICIENT !<The multiplicative coefficient applied to the RHS vector
    INTEGER(INTG), ALLOCATABLE :: RHS_DOF_TO_INTERFACE_ROW_MAP(:) !<RHS_DOF_TO_INTERFACE_ROW_MAP(rhs_dof_idx). The mapping from the rhs_dof_idx'th RHS dof in the rhs variable to the interface row.   
    INTEGER(INTG), ALLOCATABLE :: INTERFACE_ROW_TO_RHS_DOF_MAP(:) !<INTERFACE_ROW_TO_RHS_DOF_MAP(row_idx). The mapping from the row_idx'th row of the interface to the RHS dof.   
  END TYPE INTERFACE_MAPPING_RHS_TYPE
  
  TYPE INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_MATRICES !<Cache of the number of interface matrices
    INTEGER(INTG) :: LAGRANGE_VARIABLE_TYPE !<The variable type of the Lagrange field.
    REAL(DP), ALLOCATABLE :: MATRIX_COEFFICIENTS(:) !<MATRIX_COEFFICIENTS(matrix_idx). The matrix cooefficient for the matrix_idx'th interface matrix.
    LOGICAL, ALLOCATABLE :: HAS_TRANSPOSE(:) !<HAS_TRANSPOSE(matrix_idx). .TRUE. if the matrix_idx'th interface matrix has an tranpose, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: MATRIX_ROW_FIELD_VARIABLE_INDICES(:) !<MATRIX_ROW_FIELD_VARIABLE_INDICES(variable_idx). The field variable index that are mapped to the the interface matrix rows.
    INTEGER(INTG), ALLOCATABLE :: MATRIX_COL_FIELD_VARIABLE_INDICES(:) !<MATRIX_COL_FIELD_VARIABLE_INDICES(variable_idx). The field variable index that are mapped to the the interface matrix columns.
    INTEGER(INTG) :: RHS_LAGRANGE_VARIABLE_TYPE !<The Lagrange variable type mapped to the rhs vector
    REAL(DP) :: RHS_COEFFICIENT !<The coefficient multiplying the RHS vector.
  END TYPE INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE
  
  !>Contains information on an interface mapping. TODO: Generalise to non-Lagrange multipler mappings
  TYPE INTERFACE_MAPPING_TYPE  
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS  !<A pointer to the interface equations for this interface mapping
    LOGICAL :: INTERFACE_MAPPING_FINISHED !<Is .TRUE. if the interface mapping has finished being created, .FALSE. if not.
    INTEGER(INTG) :: LAGRANGE_VARIABLE_TYPE !<The variable type of the mapped Lagrange field variable.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: LAGRANGE_VARIABLE !<A pointer to the variable that is mapped to the Lagrange multiplier field variable.
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in the interface mapping
    INTEGER(INTG) :: TOTAL_NUMBER_OF_COLUMNS !<The total number of columns in the interface mapping
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_COLUMNS !<The global number of columns in the interface mapping
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOFS_MAPPING !<A pointer to the domain mapping for the columns
    INTEGER(INTG), ALLOCATABLE :: LAGRANGE_DOF_TO_COLUMN_MAP(:) !<LAGRANGE_DOF_TO_COLUMN_MAP(dof_idx). The mapping from the dof_idx'th Lagrange dof to the interface matrices column
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_MATRICES !<The number of interface matrices that the mapping is set up for.
    TYPE(INTERFACE_MATRIX_TO_VAR_MAP_TYPE), ALLOCATABLE :: INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(:) !<INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx). Information for the interface matrix rows to dependent variable maps for the interface_matrix_idx'th interface matrix.
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING !<A pointer to the interface mapping for the RHS vector.
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
  END TYPE INTERFACE_MAPPING_TYPE

  !>Contains information about the interpolation for a parameter set in interface equations
  TYPE INTERFACE_EQUATIONS_INTERPOLATION_SET_TYPE
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:) !<INTERPOLATION_PARAMETERS(field_variable_type). A pointer to the field_variable_type'th field interpolation parameters.
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:) !<INTERPOLATED_POINT(field_variable_type). A pointer to the field_variable_type'th field interpolated point. 
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: INTERPOLATED_POINT_METRICS(:) !<INTERPOLATED_POINT_METRICS(field_variable_type). A pointer to the field_variable_type'th field interpolated point metrics.
  END TYPE INTERFACE_EQUATIONS_INTERPOLATION_SET_TYPE

  !>Contains information about the interpolation for a domain (interface or coupled mesh) in the interface equations
  TYPE INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field for the domain
    INTEGER(INTG) :: NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS !<The number of geometric interpolation sets in the domain
    TYPE(INTERFACE_EQUATIONS_INTERPOLATION_SET_TYPE), ALLOCATABLE :: GEOMETRIC_INTERPOLATION(:) !<GEOMETRIC_INTERPOLATION(interpolation_set_idx). The geometric interpolation information for the interpolation_set_idx'th interpolation set.
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<A pointer to the dependent field for the domain
    INTEGER(INTG) :: NUMBER_OF_DEPENDENT_INTERPOLATION_SETS !<The number of dependent interpolation sets in the domain
    TYPE(INTERFACE_EQUATIONS_INTERPOLATION_SET_TYPE), ALLOCATABLE :: DEPENDENT_INTERPOLATION(:) !<DEPENDENT_INTERPOLATION(interpolation_set_idx). The dependent interpolation information for the interpolation_set_idx'th interpolation set.
  END TYPE INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE
  
  !>Contains information on the interpolation for the interface equations
  TYPE INTERFACE_EQUATIONS_INTERPOLATION_TYPE
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the equations
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE) :: INTERFACE_INTERPOLATION !<The interpolation information for interpolating on the interface.
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), ALLOCATABLE :: VARIABLE_INTERPOLATION(:) !<VARIABLE_INTERPOLATION(variable_idx). The interpolation for the variable_idx'th field variable in the interface condition.
  END TYPE INTERFACE_EQUATIONS_INTERPOLATION_TYPE

  !>Contains information about the interface equations for an interface condition. 
  TYPE INTERFACE_EQUATIONS_TYPE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    LOGICAL :: INTERFACE_EQUATIONS_FINISHED !<Is .TRUE. if the interface equations have finished being created, .FALSE. if not.
    INTEGER(INTG) :: OUTPUT_TYPE !<The output type for the interface equations 
    INTEGER(INTG) :: SPARSITY_TYPE !<The sparsity type for the interface equation matrices of the interface equations
    TYPE(INTERFACE_EQUATIONS_INTERPOLATION_TYPE), POINTER :: INTERPOLATION !<A pointer to the interpolation information used in the interface equations.
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface equations mapping for the interface.
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface equations matrices and vectors used for the interface equations.
  END TYPE INTERFACE_EQUATIONS_TYPE
  
  !>Contains information on the geometry for an interface condition
  TYPE INTERFACE_GEOMETRY_TYPE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition.
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<The geometric field for this equations set.
  END TYPE INTERFACE_GEOMETRY_TYPE

  !>Contains information about the Lagrange field information for an interface condition. 
  TYPE INTERFACE_LAGRANGE_TYPE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    LOGICAL :: LAGRANGE_FINISHED !<Is .TRUE. if the interface Lagrange field has finished being created, .FALSE. if not.
    LOGICAL :: LAGRANGE_FIELD_AUTO_CREATED !<Is .TRUE. if the Lagrange field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD !<A pointer to the lagrange field.
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !<The number of components in the Lagrange field.
  END TYPE INTERFACE_LAGRANGE_TYPE

  !>Contains information about the dependent field information for an interface condition. 
  TYPE INTERFACE_DEPENDENT_TYPE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    INTEGER(INTG) :: NUMBER_OF_DEPENDENT_VARIABLES !<The number of dependent variables in the interface condition
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: EQUATIONS_SETS(:) !<EQUATIONS_SETS(variable_idx). The pointer to the equations set containing the dependent variable for the variable_idx'th added dependent variable.
    TYPE(FIELD_VARIABLE_PTR_TYPE), POINTER :: FIELD_VARIABLES(:) !<FIELD_VARIABLES(variable_idx). The pointer to the variable_idx'th dependent variable in the interface condition.
    INTEGER(INTG), POINTER :: VARIABLE_MESH_INDICES(:) !<VARIABLE_MESH_INDICES(variable_idx). The mesh index of the variable_idx'th dependent variable in the interface condition.
  END TYPE INTERFACE_DEPENDENT_TYPE

  !>Contains information for the interface condition data.
  TYPE INTERFACE_CONDITION_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user identifying number of the interface condition. Must be unique
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global index of the interface condition in the interface conditions.
    LOGICAL :: INTERFACE_CONDITION_FINISHED !<Is .TRUE. ifand where  the interfaand where and where ce condition has finished being created, .FALSE. if not.
    TYPE(INTERFACE_CONDITIONS_TYPE), POINTER :: INTERFACE_CONDITIONS !<A pointer back to the interface conditions
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer back to the interface.
    INTEGER(INTG) :: METHOD !<An integer which denotes the interface condition method. \see INTERFACE_CONDITIONS_Methods,INTERFACE_CONDITIONS
    INTEGER(INTG) :: OPERATOR !<An integer which denotes whether type of interface operator. \see INTERFACE_CONDITIONS_Operator,INTERFACE_CONDITIONS
    TYPE(INTERFACE_GEOMETRY_TYPE) :: GEOMETRY !<The geometry information for the interface condition.
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE !<A pointer to the interface condition Lagrange multipler information if there are any for this interface condition.
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: DEPENDENT !<A pointer to the interface condition dependent field information if there is any for this interface condition.
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations if there are any for this interface condition
  END TYPE INTERFACE_CONDITION_TYPE

  !>A buffer type to allow for an array of pointers to a INTERFACE_CONDITION_TYPE.
  TYPE INTERFACE_CONDITION_PTR_TYPE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: PTR !<The pointer to the interface condition.
  END TYPE INTERFACE_CONDITION_PTR_TYPE
  
  !>Contains information for interface region specific data that is not of 'region' importance. <<>>
  TYPE INTERFACE_CONDITIONS_TYPE
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer back to the interface containing the interface conditions.
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_CONDITIONS !<The number of interface conditions
    TYPE(INTERFACE_CONDITION_PTR_TYPE), POINTER :: INTERFACE_CONDITIONS(:) !<INTERFACE_CONDITIONS(interface_condition_idx). A pointer to the interface_condition_idx'th interface condition.
  END TYPE INTERFACE_CONDITIONS_TYPE
  
  !>Contains information on a mesh connectivity point
  TYPE INTERFACE_ELEMENT_CONNECTIVITY_TYPE
    INTEGER(INTG) :: COUPLED_MESH_ELEMENT_NUMBER !<GLOBAL_MESH_ELEMENT_NUMBERS(connectivity_point_idx) !\todo Comment
    REAL(DP), ALLOCATABLE :: XI(:,:,:) !<XI(xi_idx,mesh_component,element_parameter_idx) !\todo Comment
  END TYPE INTERFACE_ELEMENT_CONNECTIVITY_TYPE

  !>Contains information on the coupling between meshes in an interface
  TYPE INTERFACE_MESH_CONNECTIVITY_TYPE
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer back to the interface for the coupled mesh connectivity
    TYPE(MESH_TYPE), POINTER :: INTERFACE_MESH
    TYPE(BASIS_TYPE), POINTER :: BASIS
    LOGICAL :: MESH_CONNECTIVITY_FINISHED !<Is .TRUE. if the coupled mesh connectivity has finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_INT_ELEM !<Is the number of elements within the interface mesh
    INTEGER(INTG) :: NUMBER_INT_DOM !<Is the number of domains coupled via the interface
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), ALLOCATABLE :: ELEMENTS_CONNECTIVITY(:,:) !<ELEMENTS_CONNECTIVITY(element_idx,coupled_mesh_idx)
  END TYPE INTERFACE_MESH_CONNECTIVITY_TYPE
 
  !>Contains information for the interface data.
  TYPE INTERFACE_TYPE
    INTEGER :: USER_NUMBER !<The user defined identifier for the interface. The user number must be unique.
    INTEGER :: GLOBAL_NUMBER !<The global number of the interface in the list of interfaces for a particular parent region.
    LOGICAL :: INTERFACE_FINISHED !<Is .TRUE. if the interface has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: LABEL !<A user defined label for the region.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system used by the interface. 
    TYPE(INTERFACES_TYPE), POINTER :: INTERFACES !<A pointer back to the parent interfaces
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A point to the parent region containing the interface.
    INTEGER(INTG) :: NUMBER_OF_COUPLED_MESHES !<The number of coupled meshes in the interface.
    TYPE(MESH_PTR_TYPE), POINTER :: COUPLED_MESHES(:) !<COUPLED_MESHES(mesh_idx). COUPLED_MESHES(mesh_idx)%PTR is the pointer to the mesh_idx'th mesh involved in the interface.
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: MESH_CONNECTIVITY !<A pointer to the meshes connectivity the interface.
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS  !<A pointer to the data points defined in an interface.          
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes in an interface
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the mesh in an interface.
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<A pointer to the generated meshes in an interface.
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointer to the fields defined over an interface.
    TYPE(INTERFACE_CONDITIONS_TYPE), POINTER :: INTERFACE_CONDITIONS !<The pointer to the interface conditions for this interface
  END TYPE INTERFACE_TYPE

  !>A buffer type to allow for an array of pointers to a INTERFACE_TYPE.
  TYPE INTERFACE_PTR_TYPE
    TYPE(INTERFACE_TYPE), POINTER :: PTR !<The pointer to the interface.
  END TYPE INTERFACE_PTR_TYPE
  
  !>Contains information for interfaces on a parent region.
  TYPE INTERFACES_TYPE
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer back to the parent region containing the interfaces.
    INTEGER(INTG) :: NUMBER_OF_INTERFACES !<The number of interfaces
    TYPE(INTERFACE_PTR_TYPE), POINTER :: INTERFACES(:) !<INTERFACES(interface_idx). A pointer to the interface_idx'th interface.
  END TYPE INTERFACES_TYPE

  !
  !================================================================================================================================
  !
  ! CellML types (belongs under field types?)

  !>This type is a wrapper for the C_PTR which references the actual CellML model definition object.
  TYPE CELLML_MODEL_TYPE
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment.
    INTEGER(INTG) :: GLOBAL_NUMBER !< The global number of this CellML model within the parent CellML environment.
    TYPE(VARYING_STRING) :: MODEL_ID !<The ID of the model.
    TYPE(C_PTR) :: PTR !< The handle for the actual C++ CellML model definition object
    INTEGER(INTG) :: NUMBER_OF_STATE !<The number of state variables in the CellML model.
    TYPE(VARYING_STRING), ALLOCATABLE :: STATE_VARIABLE_ID(:) !<STATE_VARIABLE_ID(state_variable_idx). The ID for the state_variable_idx'th state variable.
    INTEGER(INTG) :: NUMBER_OF_INTERMEDIATE !<The number of intermediate variables in the CellML model.
    TYPE(VARYING_STRING), ALLOCATABLE :: INTERMEDIATE_VARIABLE_ID(:) !<INTERMEDIATE_VARIABLE_ID(intermediate_variable_idx). The ID for the intermediate_variable_idx'th intermediate variable.
    INTEGER(INTG) :: NUMBER_OF_PARAMETERS !<The number of parameters in the CellML model.
    TYPE(VARYING_STRING), ALLOCATABLE :: PARAMETER_VARIABLE_ID(:) !<PARAMETER_VARIABLE_ID(parameter_variable_idx). The ID for the parameter_variable_idx'th parameter variable.
 END TYPE CELLML_MODEL_TYPE

  !>A buffer type to allow for an array of pointers to a CELLML_MODEL_TYPE
  TYPE CELLML_MODEL_PTR_TYPE
    TYPE(CELLML_MODEL_TYPE), POINTER :: PTR
  END TYPE CELLML_MODEL_PTR_TYPE

  !>Contains information on the models field for a CellML environment
  TYPE CELLML_MODELS_FIELD_TYPE
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment 
    LOGICAL :: MODELS_FIELD_FINISHED  !<Is .TRUE. if the models field has finished being created, .FALSE. if not.
    LOGICAL :: MODELS_FIELD_AUTO_CREATED !<Is .TRUE. if the models field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: MODELS_FIELD !<A pointer to the models field
    INTEGER(INTG) :: ONLY_ONE_MODEL_INDEX !<If only one model is used in the models field for the CellML environment then this will be equal to the model index. It will be zero otherwise.
  END TYPE CELLML_MODELS_FIELD_TYPE
  
  !>Contains information on the state field for a CellML environment
  TYPE CELLML_STATE_FIELD_TYPE
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment 
    LOGICAL :: STATE_FIELD_FINISHED  !<Is .TRUE. if the state field has finished being created, .FALSE. if not.
    LOGICAL :: STATE_FIELD_AUTO_CREATED !<Is .TRUE. if the state field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: STATE_FIELD !<A pointer to the state field
  END TYPE CELLML_STATE_FIELD_TYPE
  
  !>Contains information on the intermediate field for a CellML environment
  TYPE CELLML_INTERMEDIATE_FIELD_TYPE
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment 
    LOGICAL :: INTERMEDIATE_FIELD_FINISHED  !<Is .TRUE. if the intermediate field has finished being created, .FALSE. if not.
    LOGICAL :: INTERMEDIATE_FIELD_AUTO_CREATED !<Is .TRUE. if the intermediate field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: INTERMEDIATE_FIELD !<A pointer to the intermediate field
  END TYPE CELLML_INTERMEDIATE_FIELD_TYPE
  
  !>Contains information on the parameters field for a CellML environment
  TYPE CELLML_PARAMETERS_FIELD_TYPE
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment 
    LOGICAL :: PARAMETERS_FIELD_FINISHED  !<Is .TRUE. if the parameters field has finished being created, .FALSE. if not.
    LOGICAL :: PARAMETERS_FIELD_AUTO_CREATED !<Is .TRUE. if the parameters field has been auto created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: PARAMETERS_FIELD !<A pointer to the parameters field
  END TYPE CELLML_PARAMETERS_FIELD_TYPE
  
  !>Contains information on the mapping between CellML fields and OpenCMISS fields and vise versa.
  TYPE CELLML_MODEL_MAP_TYPE
    INTEGER(INTG) :: CELLML_MAP_TYPE !<The direction of the mapping. \see CELLML_FieldMappingTypes,CMISS_CELLML
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to/from the field being mapped.
    INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type being mapped.
    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number being mapped.
    INTEGER(INTG) :: FIELD_PARAMETER_SET !<The field variable parameter set being mapped.
    TYPE(VARYING_STRING) :: VARIABLE_ID !<The variable ID of the CellML variable being mapped.
    INTEGER(INTG) :: CELLML_FIELD_TYPE !<The type of CellML field variable being mapped. \see CELLML_FieldTypes,CMISS_CELLML
    INTEGER(INTG) :: CELLML_VARIABLE_NUMBER !<The CellML variable component number being mapped.
    INTEGER(INTG) :: CELLML_PARAMETER_SET !<The CellML variable parameter set being mapped.
  END TYPE CELLML_MODEL_MAP_TYPE

  !>Buffer type to allow an array of pointers to CELLML_MODEL_MAP_FIELD_TYPE
  TYPE CELLML_MODEL_MAP_PTR_TYPE
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: PTR !<A pointer to the CELLML_MODEL_MAP_TYPE
  END TYPE CELLML_MODEL_MAP_PTR_TYPE

  !>Contains information on the maps between a CellML model and external OpenCMISS fields.
  TYPE CELLML_MODEL_MAPS_TYPE
    INTEGER(INTG) :: NUMBER_OF_FIELDS_MAPPED_TO !<The number of OpenCMISS fields mapped to this CellML model's variables.
    TYPE(CELLML_MODEL_MAP_PTR_TYPE), ALLOCATABLE :: FIELDS_MAPPED_TO(:) !<FIELDS_MAPPED_TO(map_idx). The map_idx'th field mapping for OpenCMISS fields that the CellML model maps to.
    INTEGER(INTG) :: NUMBER_OF_FIELDS_MAPPED_FROM !<The number of CellML variable fields mapped to OpenCMISS fields.
    TYPE(CELLML_MODEL_MAP_PTR_TYPE), ALLOCATABLE :: FIELDS_MAPPED_FROM(:) !<FIELDS_MAPPED_FROM(map_idx). The map_idx'th field mapping for OpenCMISS fields that the CellML model maps from.
  END TYPE CELLML_MODEL_MAPS_TYPE

  !>Buffer type to allow arrays of pointer to CELLML_MODEL_MAPS_TYPE
  TYPE CELLML_MODEL_MAPS_PTR_TYPE
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: PTR !<A pointer to the CELLML_MODEL_MAPS_TYPE
  END TYPE CELLML_MODEL_MAPS_PTR_TYPE
  
  !>Contains information on the maps between CellML and external OpenCMISS fields.
  TYPE CELLML_FIELD_MAPS_TYPE
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer back to the CellML environment.
    LOGICAL :: CELLML_FIELD_MAPS_FINISHED !<Is .TRUE. if the CellML maps have finished being created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: SOURCE_GEOMETRIC_FIELD !<The source geometric field for the CellML environment.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: SOURCE_FIELD_VARIABLE !<The source field variable for the CellML environment.
    TYPE(DOMAIN_TYPE), POINTER :: SOURCE_FIELD_DOMAIN !<The source field domain for the CellML environment.
    INTEGER(INTG) :: SOURCE_FIELD_INTERPOLATION_TYPE !<The source field interpolation type for the CellML environment.
    TYPE(CELLML_MODEL_MAPS_PTR_TYPE), ALLOCATABLE :: MODEL_MAPS(:) !<MODEL_MAPS(model_idx). Contains information on the maps between the model_idx'th CellML model and external OpenCMISS fields.
    !INTEGER(INTG) :: NUMBER_OF_SOURCE_DOFS !<The number of local (excluding ghosts) source dofs.
    !INTEGER(INTG) :: TOTAL_NUMBER_OF_SOURCE_DOFS !<The number of local (including ghosts) source dofs.
    !INTEGER(INTG) :: GLOBAL_NUMBER_OF_SOURCE_DOFS !<The number of global source dofs.
  END TYPE CELLML_FIELD_MAPS_TYPE
  
  !>Contains information for a CellML environment.
  TYPE CELLML_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing this CellML environment.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the CellML environment in the list of environments for a field.
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the CellML environment. The user number must be unique.
    TYPE(CELLML_ENVIRONMENTS_TYPE), POINTER :: ENVIRONMENTS !<A pointer back to the CellML environments.
    LOGICAL :: CELLML_FINISHED !<Is .TRUE. if the environment has finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_MODELS !< The number of models defined in the CellML environment
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_STATE !<The maximum number of state variables across all models.
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_PARAMETERS !<The maximum number of parameters variables across all models.
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_INTERMEDIATE !<The maximum number of intermediate variables across all models.
    TYPE(CELLML_MODEL_PTR_TYPE), ALLOCATABLE :: MODELS(:) !< MODELS(model_idx). The array of pointers to the models.
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: FIELD_MAPS !<A pointer the information on CellML<-->Field maps.
    TYPE(CELLML_MODELS_FIELD_TYPE), POINTER :: MODELS_FIELD !<A pointer to the models field information.
    TYPE(CELLML_STATE_FIELD_TYPE), POINTER :: STATE_FIELD !<A pointer to the state field information.
    TYPE(CELLML_INTERMEDIATE_FIELD_TYPE), POINTER :: INTERMEDIATE_FIELD !<A pointer to the intermediate field information.
    TYPE(CELLML_PARAMETERS_FIELD_TYPE), POINTER :: PARAMETERS_FIELD !<A pointer to the parameters field information.
    LOGICAL :: CELLML_GENERATED !<Is .TRUE. if the CellML environment has finished being generated, .FALSE. if not.
  END TYPE CELLML_TYPE

  !> A buffer type to allow for an array of pointers to a CELLML_TYPE.
  !! \todo Is this needed? not currently used...
  TYPE CELLML_PTR_TYPE
    TYPE(CELLML_TYPE), POINTER :: PTR !< The pointer to the CellML environment.
  END TYPE CELLML_PTR_TYPE

  !>Contains information on the CellML environments defined.
  TYPE CELLML_ENVIRONMENTS_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the CellML environments
    INTEGER(INTG) :: NUMBER_OF_ENVIRONMENTS !<The number of environments defined.
    TYPE(CELLML_PTR_TYPE), ALLOCATABLE :: ENVIRONMENTS(:) !<The array of pointers to the CellML environments.
  END TYPE CELLML_ENVIRONMENTS_TYPE

  !
  !================================================================================================================================
  !
  ! Solver matrices types
  !
  
  !>Contains information on the solver matrix
  TYPE SOLVER_MATRIX_TYPE
    INTEGER(INTG) :: MATRIX_NUMBER !<The number of the solver matrix
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices for this solver matrix
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if the solver matrix is to be updated
    INTEGER(INTG) :: STORAGE_TYPE !<The storage type for the solver matrix.
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in the distributed solver matrix
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOLVER_VECTOR !<A pointer to the distributed solver vector associated with the matrix
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the distributed solver matrix data
  END TYPE SOLVER_MATRIX_TYPE

  !>A buffer type to allow for an array of pointers to a SOLVER_MATRIX_TYPE \see TYPES:SOLUTION_MATRIX_TYPE
  TYPE SOLVER_MATRIX_PTR_TYPE
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: PTR !<The pointer to the solver matrix.
  END TYPE SOLVER_MATRIX_PTR_TYPE
  
  !>Contains information on the solver matrices and rhs vector
  TYPE SOLVER_MATRICES_TYPE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver
    LOGICAL :: SOLVER_MATRICES_FINISHED !<Is .TRUE. if the solver matrices have finished being created, .FALSE. if not.
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping for these solver matrices
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of (local) rows in the distributed solution matrix for this computational node
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ROWS !<The number of global rows in the distributed solution matrix
    INTEGER(INTG) :: LIBRARY_TYPE !<The library type for the solver matrices
    !Linear matrices
    INTEGER(INTG) :: NUMBER_OF_MATRICES !<The number of solver matrices defined for the problem
    TYPE(SOLVER_MATRIX_PTR_TYPE), ALLOCATABLE :: MATRICES(:) !<MATRICES(matrix_idx)%PTR contains the information on the matrix_idx'th solver matrix
    !Nonlinear matrices and vectors
    LOGICAL :: UPDATE_RESIDUAL !<Is .TRUE. if the residual vector is to be updated
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL !<A pointer to the distributed residual vector for nonlinear problems
    !Right hand side vector
    LOGICAL :: UPDATE_RHS_VECTOR !<Is .TRUE. if the RHS vector is to be updated
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RHS_VECTOR !<A pointer to the distributed RHS vector for the solver matrices
  END TYPE SOLVER_MATRICES_TYPE

  !
  !================================================================================================================================
  !
  ! Solver equations types
  !

  !>Contains information about the solver equations for a solver. \see OPENCMISS::CMISSSolverEquationsType
  TYPE SOLVER_EQUATIONS_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    LOGICAL :: SOLVER_EQUATIONS_FINISHED !<Is .TRUE. if the solver equations have finished being created, .FALSE. if not.

    INTEGER(INTG) :: LINEARITY !<The linearity type of the solver equations
    INTEGER(INTG) :: TIME_DEPENDENCE !<The time dependence type of the solver equations

    INTEGER(INTG) :: SPARSITY_TYPE !<The type of sparsity to use in the solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES

    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices for the problem

    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary condition information for the solver equations.

  END TYPE SOLVER_EQUATIONS_TYPE

  !
  !================================================================================================================================
  !
  ! CellML equations types
  !

  !>Contains information about the CellML equations for a solver.
  TYPE CELLML_EQUATIONS_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    LOGICAL :: CELLML_EQUATIONS_FINISHED !<Is .TRUE. if the CellML equations have finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_CELLML_ENVIRONMENTS !<The number of CellML environments in the equations
    TYPE(CELLML_PTR_TYPE), ALLOCATABLE :: CELLML_ENVIRONMENTS(:) !<CELLML_ENVIORNMENTS(cellml_idx). The array of pointers to the CellML environments for these CellML equations.     
  END TYPE CELLML_EQUATIONS_TYPE
  
  !
  !================================================================================================================================
  !
  ! Solver types
  !
  
  !>Contains information for a dynamic solver
  TYPE DYNAMIC_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the dynamic solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    LOGICAL :: SOLVER_INITIALISED !<Is .TRUE. if the solver has been initialised, .FALSE. if not.
    INTEGER(INTG) :: LINEARITY !<The linearity type of the dynamic solver \see SOLVER_ROUTINES_DynamicLinearityTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: ORDER !<The order of the dynamic solve \see SOLVER_ROUTINES_DynamicOrderTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: DEGREE !<The degree of the time interpolation polynomial \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: SCHEME !<The dyanamic solver scheme \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
    REAL(DP), ALLOCATABLE :: THETA(:) !<THETA(degree_idx). The theta value for the degree_idx'th polynomial in the dynamic solver
    LOGICAL :: EXPLICIT !<Is .TRUE. if the dynamic scheme is an explicit scheme, .FALSE. if not.
    LOGICAL :: RESTART !<Is .TRUE. if the dynamic scheme is to be restarted (i.e., recalculate values at the current time step), .FALSE. if not.
    LOGICAL :: ALE !<Is .TRUE. if the dynamic scheme is an ALE scheme, .FALSE. if not.
    LOGICAL :: UPDATE_BC !<Is .TRUE. if the dynamic scheme has changing bc, .FALSE. if not.
    REAL(DP) :: CURRENT_TIME !<The current time value for the dynamic solver.
    REAL(DP) :: TIME_INCREMENT !<The time increment for the dynamic solver to solver for.
    TYPE(SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linked linear solver
    TYPE(SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer to the linked nonlinear solver
  END TYPE DYNAMIC_SOLVER_TYPE
  
   !>Contains information for an forward Euler differential-algebraic equation solver
  TYPE FORWARD_EULER_DAE_SOLVER_TYPE
    TYPE(EULER_DAE_SOLVER_TYPE), POINTER :: EULER_DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the forward Euler differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE FORWARD_EULER_DAE_SOLVER_TYPE

  !>Contains information for an backward Euler differential-algebraic equation solver
  TYPE BACKWARD_EULER_DAE_SOLVER_TYPE
    TYPE(EULER_DAE_SOLVER_TYPE), POINTER :: EULER_DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the backward Euler differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE BACKWARD_EULER_DAE_SOLVER_TYPE

  !>Contains information for an improved Euler differential-algebraic equation solver
  TYPE IMPROVED_EULER_DAE_SOLVER_TYPE
    TYPE(EULER_DAE_SOLVER_TYPE), POINTER :: EULER_DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the improved Euler differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE IMPROVED_EULER_DAE_SOLVER_TYPE
  
  !>Contains information for an Euler differential-algebraic equation solver
  TYPE EULER_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: EULER_TYPE !<The type of Euler differential-algebraic equation solver \see SOLVER_ROUTINES_EulerDESolverTypes,SOLVER_ROUTINES
    TYPE(FORWARD_EULER_DAE_SOLVER_TYPE), POINTER :: FORWARD_EULER_SOLVER !<A pointer to the forward Euler solver information
    TYPE(BACKWARD_EULER_DAE_SOLVER_TYPE), POINTER :: BACKWARD_EULER_SOLVER !<A pointer to the backward Euler solver information
    TYPE(IMPROVED_EULER_DAE_SOLVER_TYPE), POINTER :: IMPROVED_EULER_SOLVER !<A pointer to the improved Euler solver information
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the Euler differential equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE EULER_DAE_SOLVER_TYPE

  !>Contains information for a Crank-Nicholson differential-algebraic equation solver
  TYPE CRANK_NICOLSON_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the Crank-Nicholson differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE CRANK_NICOLSON_DAE_SOLVER_TYPE
  
  !>Contains information for a Runge-Kutta differential-algebraic equation solver
  TYPE RUNGE_KUTTA_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the Runge-Kutta differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE RUNGE_KUTTA_DAE_SOLVER_TYPE
  
  !>Contains information for an Adams-Moulton differential-algebraic equation solver
  TYPE ADAMS_MOULTON_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the Adams-Moulton differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE ADAMS_MOULTON_DAE_SOLVER_TYPE
  
  !>Contains information for a BDF differential-algebraic equation solver
  TYPE BDF_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the BDF differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE BDF_DAE_SOLVER_TYPE
  
  !>Contains information for a Rush-Larson differential-algebraic equation solver
  TYPE RUSH_LARSON_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the Rush-Larson differential-algebraic equation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  END TYPE RUSH_LARSON_DAE_SOLVER_TYPE
  
  !>Contains information for an external differential-algebraic equation solver
  TYPE EXTERNAL_DAE_SOLVER_TYPE
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic solver
  END TYPE EXTERNAL_DAE_SOLVER_TYPE
  
  !>Contains information for an differential-algebraic equation solver
  TYPE DAE_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG) :: DAE_TYPE !<The differential-algebraic equation type \see SOLVER_ROUTINES_DAETypes,SOLVER_ROUTINES
    INTEGER(INTG) :: DAE_SOLVE_TYPE !<The solve type for the differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolveTypes,SOLVER_ROUTINES
    REAL(DP) :: START_TIME !<The start time to integrate from
    REAL(DP) :: END_TIME !<The end time to integrate to
    REAL(DP) :: INITIAL_STEP !<The (initial) time step
    TYPE(EULER_DAE_SOLVER_TYPE), POINTER :: EULER_SOLVER !<A pointer to information for an Euler solver
    TYPE(CRANK_NICOLSON_DAE_SOLVER_TYPE), POINTER :: CRANK_NICOLSON_SOLVER !<A pointer to information for a Crank-Nicholson solver
    TYPE(RUNGE_KUTTA_DAE_SOLVER_TYPE), POINTER :: RUNGE_KUTTA_SOLVER !<A pointer to information for a Runge-Kutta solver
    TYPE(ADAMS_MOULTON_DAE_SOLVER_TYPE), POINTER :: ADAMS_MOULTON_SOLVER !<A pointer to information for an Adams-Moulton solver
    TYPE(BDF_DAE_SOLVER_TYPE), POINTER :: BDF_SOLVER !<A pointer to information for a BDF solver
    TYPE(RUSH_LARSON_DAE_SOLVER_TYPE), POINTER :: RUSH_LARSON_SOLVER !<A pointer to information for a Rush-Larson solver
    TYPE(EXTERNAL_DAE_SOLVER_TYPE), POINTER :: EXTERNAL_SOLVER !<A pointer to information for an external solver
  END TYPE DAE_SOLVER_TYPE
  
  !>Contains information for a direct linear solver
  TYPE LINEAR_DIRECT_SOLVER_TYPE
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linear solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the linear direct solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLVER_MATRICES_LIBRARY !<The library type for the linear direct solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG) :: DIRECT_SOLVER_TYPE !<The type of direct linear solver
    TYPE(PETSC_PC_TYPE) :: PC !<The PETSc preconditioner object
    TYPE(PETSC_KSP_TYPE) :: KSP !<The PETSc solver object
  END TYPE LINEAR_DIRECT_SOLVER_TYPE

  !>Contains information for an iterative linear solver
  TYPE LINEAR_ITERATIVE_SOLVER_TYPE
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linear solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the linear solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLVER_MATRICES_LIBRARY !<The library type for the linear iterative solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG) :: ITERATIVE_SOLVER_TYPE !<The type of iterative solver \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: ITERATIVE_PRECONDITIONER_TYPE !<The type of iterative preconditioner \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLUTION_INITIALISE_TYPE !<The type of solution vector initialisation \see SOLVER_ROUTINES_SolutionInitialiseTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ITERATIONS !<The maximum number of iterations
    REAL(DP) :: RELATIVE_TOLERANCE !<The relative tolerance between the rhs and residual norm
    REAL(DP) :: ABSOLUTE_TOLERANCE !<The absolute tolerance of the residual norm
    REAL(DP) :: DIVERGENCE_TOLERANCE !<The absolute tolerance of the residual norm
    INTEGER(INTG) :: GMRES_RESTART !<The GMRES restart iterations size
    TYPE(PETSC_PC_TYPE) :: PC !<The PETSc preconditioner object
    TYPE(PETSC_KSP_TYPE) :: KSP !<The PETSc solver object
  END TYPE LINEAR_ITERATIVE_SOLVER_TYPE
  
  !>Contains information for a linear solver
  TYPE LINEAR_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG) :: LINEAR_SOLVE_TYPE !<The type of linear solver \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
    LOGICAL :: LINKED_NEWTON_PETSC_SOLVER !<Is .TRUE. if this linear solver has been linked from a PETSc Newton solver, .FALSE. if not.
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: DIRECT_SOLVER !<A pointer to the direct solver information
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: ITERATIVE_SOLVER !<A pointer to the iterative solver information
  END TYPE LINEAR_SOLVER_TYPE

  !>Contains information for a Newton line search nonlinear solver
  TYPE NEWTON_LINESEARCH_SOLVER_TYPE
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer to the Newton solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the Newton linesearch solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLVER_MATRICES_LIBRARY !<The library type for the Newton linesearch solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG) :: LINESEARCH_TYPE !<The line search type \see SOLVER_ROUTINES_NonlinearLineSearchTypes,SOLVER_ROUTINES
    REAL(DP) :: LINESEARCH_ALPHA !<The line search alpha
    REAL(DP) :: LINESEARCH_MAXSTEP !<The line search maximum step
    REAL(DP) :: LINESEARCH_STEPTOLERANCE !<The line search step tolerance
    TYPE(PETSC_ISCOLORING_TYPE) :: JACOBIAN_ISCOLORING !<The Jacobian matrix index set colouring
    TYPE(PETSC_MATFDCOLORING_TYPE) :: JACOBIAN_FDCOLORING !<The Jacobian matrix finite difference colouring
    TYPE(PETSC_SNES_TYPE) :: SNES !<The PETSc nonlinear solver object
  END TYPE NEWTON_LINESEARCH_SOLVER_TYPE
  
  !>Contains information for a Newton trust region nonlinear solver
  TYPE NEWTON_TRUSTREGION_SOLVER_TYPE
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer to the Newton solver 
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the nonlinear solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLVER_MATRICES_LIBRARY !<The library type for the Newton trustregion solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    REAL(DP) :: TRUSTREGION_TOLERANCE !<The trust region tolerance
    REAL(DP) :: TRUSTREGION_DELTA0 !<The trust region delta0
    TYPE(PETSC_SNES_TYPE) :: SNES !<The PETSc nonlinear solver object
  END TYPE NEWTON_TRUSTREGION_SOLVER_TYPE

  !>Contains information for a Newton nonlinear solver
  TYPE NEWTON_SOLVER_TYPE
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer to the nonlinear solver
    INTEGER(INTG) :: NEWTON_SOLVE_TYPE !<The type of Newton solver
    INTEGER(INTG) :: SOLUTION_INITIALISE_TYPE !<The type of solution vector initialisation \see SOLVER_ROUTINES_SolutionInitialiseTypes,SOLVER_ROUTINES
    INTEGER(INTG) :: TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS !<The number of function evaluations performed by the nonlinear solver
    INTEGER(INTG) :: TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS !<The number of Jacobian evaluations performed by the nonlinear solver
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ITERATIONS !<The maximum number of iterations
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_FUNCTION_EVALUATIONS !<The maximum number of function evaluations
    INTEGER(INTG) :: JACOBIAN_CALCULATION_TYPE !<The type of calculation used for the Jacobian \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
    REAL(DP) :: ABSOLUTE_TOLERANCE !<The tolerance between the absolute decrease between the solution norm and the initial guess
    REAL(DP) :: RELATIVE_TOLERANCE !<The tolerance between the relative decrease between the solution norm and the initial guess
    REAL(DP) :: SOLUTION_TOLERANCE !<The tolerance of the change in the norm of the solution
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER !<A pointer to the Newton line search solver information
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER !<A pointer to the Newton trust region solver information
    TYPE(SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linked linear solver
    TYPE(SOLVER_TYPE), POINTER :: CELLML_EVALUATOR_SOLVER !<A pointer to the linked CellML solver
  END TYPE NEWTON_SOLVER_TYPE
    
  !>Contains information for a nonlinear solver
  TYPE NONLINEAR_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the problem_solver
    INTEGER(INTG) :: NONLINEAR_SOLVE_TYPE !<The type of nonlinear solver \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer to the Newton solver information
  END TYPE NONLINEAR_SOLVER_TYPE
  
  !>Contains information for an eigenproblem solver
  TYPE EIGENPROBLEM_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the eigenproblem solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLVER_MATRICES_LIBRARY !<The library type for the eigenproblem solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
  END TYPE EIGENPROBLEM_SOLVER_TYPE
  
  !>Contains information for an optimiser solver
  TYPE OPTIMISER_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the optimiser solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG) :: SOLVER_MATRICES_LIBRARY !<The library type for the optimiser solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
  END TYPE OPTIMISER_SOLVER_TYPE

  !>Contains information for a CellML evaluation solver
  TYPE CELLML_EVALUATOR_SOLVER_TYPE
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG) :: SOLVER_LIBRARY !<The library type for the CellML evaluation solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment for the solver
    REAL(DP) :: CURRENT_TIME !<The current time value for the evaluator solver
  END TYPE CELLML_EVALUATOR_SOLVER_TYPE

   !>A buffer type to allow for an array of pointers to a SOLVER_TYPE \see TYPES::SOLVER_TYPE
  TYPE SOLVER_PTR_TYPE
    TYPE(SOLVER_TYPE), POINTER :: PTR
  END TYPE SOLVER_PTR_TYPE

 !>Contains information on the type of solver to be used. \see OPENCMISS::CMISSSolverType
  TYPE SOLVER_TYPE
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the control loop solvers
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the solver in the list of solvers
    TYPE(SOLVER_TYPE), POINTER :: LINKING_SOLVER !<A pointer to any solver that is linking to this solver
    INTEGER(INTG) :: NUMBER_OF_LINKED_SOLVERS !<The number of linked solvers
    TYPE(SOLVER_PTR_TYPE), ALLOCATABLE :: LINKED_SOLVERS(:) !<LINKED_SOLVERS(linked_solver_idx). A pointer to the linked solvers
    TYPE(SOLVER_PTR_TYPE), ALLOCATABLE :: LINKED_SOLVER_TYPE_MAP(:) !<LINKED_SOLVER_TYPE_MAP(linked_solver_idx). The map from the available linked solver types to the linked solver types that are defined for the solver. linked_solver_idx varies from 1 to SOLVER_ROUTINES::SOLVER_NUMBER_OF_SOLVER_TYPES. If the particular linked solver type has not been defined on the solver then the LINKED_SOLVER_TYPE_MAP will be NULL. \see SOLVER_ROUTINES_SolverTypes
    LOGICAL :: SOLVER_FINISHED !<Is .TRUE. if the solver has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: LABEL !<A user defined label for the solver.
     
    INTEGER(INTG) :: OUTPUT_TYPE !<The type of output required \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
    
    INTEGER(INTG) :: SOLVE_TYPE !<The type of the solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
    
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linear solver information
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer to the nonlinear solver information
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER !<A pointer to the dynamic solver information
    TYPE(DAE_SOLVER_TYPE), POINTER :: DAE_SOLVER !<A pointer to the differential-algebraic equation solver information
    TYPE(EIGENPROBLEM_SOLVER_TYPE), POINTER :: EIGENPROBLEM_SOLVER !<A pointer to the eigenproblem solver information
    TYPE(OPTIMISER_SOLVER_TYPE), POINTER :: OPTIMISER_SOLVER !<A pointer to the optimiser solver information
    TYPE(CELLML_EVALUATOR_SOLVER_TYPE), POINTER :: CELLML_EVALUATOR_SOLVER !<A pointer to the CellML solver information

    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS !<A pointer to the CellML equations
    
  END TYPE SOLVER_TYPE

  !>Contains information on the solvers to be used in a control loop
  TYPE SOLVERS_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the parent control loop
    LOGICAL :: SOLVERS_FINISHED !<Is .TRUE. if the solvers have finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_SOLVERS !<The number of solvers
    TYPE(SOLVER_PTR_TYPE), ALLOCATABLE :: SOLVERS(:) !<A pointer to the solvers information for the problem.
  END TYPE SOLVERS_TYPE
  
  !
  !================================================================================================================================
  !
  ! Solver mapping types
  !
  
  TYPE EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_SOLVER_COLS !<The number of solver cols this equations column is mapped to
    INTEGER(INTG), ALLOCATABLE :: SOLVER_COLS(:) !<SOLVER_COLS(i). The solver column number for the i'th solver column that this column is mapped to
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The coupling coefficients for the i'th solver column that this column is mapped to.
  END TYPE EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE

  TYPE EQUATIONS_TO_SOLVER_MAPS_TYPE
    INTEGER(INTG) :: EQUATIONS_MATRIX_TYPE !<The type of equations matrix \see SOLVER_MAPPING_EquationsMatrixTypes,SOLVER_MAPPING
    INTEGER(INTG) :: EQUATIONS_MATRIX_NUMBER !<The equations matrix number being mapped.
    INTEGER(INTG) :: SOLVER_MATRIX_NUMBER !<The solver matrix number being mapped.
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX !<A pointer to the equations matrix being mapped.
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix being mapped.
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE), ALLOCATABLE :: EQUATIONS_COL_TO_SOLVER_COLS_MAP(:) !<EQUATIONS_COL_TO_SOLVER_COL_MAP(column_idx). The mapping from the column_idx'th column of this equations matrix to the solver matrix columns. Note for interface matrices this will actually be for the transposed interface matrix i.e., for the interface matrix rows.
  END TYPE EQUATIONS_TO_SOLVER_MAPS_TYPE

  !>A buffer type to allow for an array of pointers to a EQUATIONS_TO_SOLVER_MAPS_TYPE \see TYPES::EQUATIONS_TO_SOLVER_MAPS_TYPE
  TYPE EQUATIONS_TO_SOLVER_MAPS_PTR_TYPE
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: PTR !<A pointer to the equations to solver maps
  END TYPE EQUATIONS_TO_SOLVER_MAPS_PTR_TYPE
  
  TYPE INTERFACE_TO_SOLVER_MAPS_TYPE
    INTEGER(INTG) :: INTERFACE_MATRIX_NUMBER !<The interface matrix number being mapped.
    INTEGER(INTG) :: SOLVER_MATRIX_NUMBER !<The solver matrix number being mapped.
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix being mapped.
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix being mapped.
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE), ALLOCATABLE :: INTERFACE_ROW_TO_SOLVER_COLS_MAP(:) !<INTERFACE_ROW_TO_SOLVER_COL_MAP(column_idx). The mapping from the row_idx'th row of this interface matrix to the solver matrix columns.
  END TYPE INTERFACE_TO_SOLVER_MAPS_TYPE

  !>A buffer type to allow for an array of pointers to a INTERFACE_TO_SOLVER_MAPS_TYPE \see TYPES::INTERFACE_TO_SOLVER_MAPS_TYPE
  TYPE INTERFACE_TO_SOLVER_MAPS_PTR_TYPE
    TYPE(INTERFACE_TO_SOLVER_MAPS_TYPE), POINTER :: PTR !<A pointer to the interface to solver maps
  END TYPE INTERFACE_TO_SOLVER_MAPS_PTR_TYPE
  
  TYPE JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_SOLVER_COLS !<The number of solver cols this Jacobian column is mapped to
    INTEGER(INTG), ALLOCATABLE :: SOLVER_COLS(:) !<SOLVER_COLS(i). The solver column number for the i'th solver column that this column is mapped to
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The coupling coefficients for the i'th solver column that this column is mapped to.
  END TYPE JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE

  TYPE JACOBIAN_TO_SOLVER_MAP_TYPE
    INTEGER(INTG) :: SOLVER_MATRIX_NUMBER !<The solver matrix number being mapped.
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX !<A pointer to the Jacobian matrix being mapped.
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix being mapped.
    TYPE(JACOBIAN_COL_TO_SOLVER_COLS_MAP_TYPE), ALLOCATABLE :: JACOBIAN_COL_TO_SOLVER_COLS_MAP(:) !<JACOBIAN_COL_TO_SOLVER_COL_MAP(column_idx). The mapping from the column_idx'th column of the Jacobian matrix to the solver matrix columns.
  END TYPE JACOBIAN_TO_SOLVER_MAP_TYPE

  TYPE JACOBIAN_TO_SOLVER_MAP_PTR_TYPE
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: PTR !<A pointer to the Jacobian to solver map
  END TYPE JACOBIAN_TO_SOLVER_MAP_PTR_TYPE

  !>Contains information on the mappings between field variable dofs inequations and the solver matrix columns (solver dofs) \todo rename solver col to be solver dof here???
  TYPE VARIABLE_TO_SOLVER_COL_MAP_TYPE
    INTEGER(INTG), ALLOCATABLE :: COLUMN_NUMBERS(:) !<COLUMN_NUMBERS(variable_dof_idx). The solver column number (solver dof) that the variable_dof_idx'th variable dof is mapped to.
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(variable_dof_idx). The multiplicative constant for the mapping between the variable_dof_idx'th variable dof and the solver dof.
    REAL(DP), ALLOCATABLE :: ADDITIVE_CONSTANTS(:) !<ADDITIVE_CONSTANTS(variable_dof_idx). The additive constant for the mapping between the variable_dof_idx'th variable dof and the solver dof.
  END TYPE VARIABLE_TO_SOLVER_COL_MAP_TYPE

  TYPE EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE_TYPE
    INTEGER(INTG) :: INTERFACE_CONDITION_INDEX
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    INTEGER(INTG) :: INTERFACE_MATRIX_NUMBER
  END TYPE EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE_TYPE
  
  !>Contains information on the equations to solver matrix mappings when indexing by solver matrix number
  TYPE EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE
    INTEGER(INTG) :: SOLVER_MATRIX_NUMBER !<The number of the solver matrix for these mappings

    INTEGER(INTG) :: NUMBER_OF_VARIABLES !<The number of variables mapped to this solver matrix
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). The variable type for the variable_idx'th variable that is mapped to the solver matrix.
    TYPE(FIELD_VARIABLE_PTR_TYPE), ALLOCATABLE :: VARIABLES(:) !<VARIABLES(variable_idx). VARIABLES(variable_idx)%PTR points to the variable_idx'th variable that is mapped to the solver matrix.
    TYPE(VARIABLE_TO_SOLVER_COL_MAP_TYPE), ALLOCATABLE :: VARIABLE_TO_SOLVER_COL_MAPS(:) !<VARIABLE_TO_SOLVER_COL_MAPS(variable_idx). The mappings from the variable dofs to the solver dofs for the variable_idx'th variable in the equations set that is mapped to the solver matrix.
    INTEGER(INTG) :: NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES !<The number of dynamic equations matrices mapped to this solver matrix
    TYPE(EQUATIONS_TO_SOLVER_MAPS_PTR_TYPE), ALLOCATABLE :: DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(:) !<DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx). The maps from the equations_idx'th linear equations matrix to solver matrix
    INTEGER(INTG) :: NUMBER_OF_LINEAR_EQUATIONS_MATRICES !<The number of linear equations matrices mapped to this solver matrix
    TYPE(EQUATIONS_TO_SOLVER_MAPS_PTR_TYPE), ALLOCATABLE :: LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(:) !<LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx). The maps from the equations_idx'th linear equations matrix to solver matrix
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS_JACOBIANS !<The number of nonlinear equations Jacobian matrices mapped to this solver matrix
    TYPE(JACOBIAN_TO_SOLVER_MAP_PTR_TYPE), ALLOCATABLE :: JACOBIAN_TO_SOLVER_MATRIX_MAPS(:) !<JACOBIAN_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx). The map from the equations_matrix_idx'th Jacobian matrix to the solver matrix
  END TYPE EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE

  !>Contains information on the equations to solver matrix mappings when indexing by equations matrix number
  TYPE EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE
    INTEGER(INTG) :: EQUATIONS_MATRIX_NUMBER !<The number of the equations matrix for these mappings
    INTEGER(INTG) :: NUMBER_OF_SOLVER_MATRICES !<The number of solver matrices mapped to this equations matrix
    TYPE(EQUATIONS_TO_SOLVER_MAPS_PTR_TYPE), ALLOCATABLE :: EQUATIONS_TO_SOLVER_MATRIX_MAPS(:) !<EQUATIONS_TO_SOLVER_MATRIX_MAPS(solver_matrix_idx). The maps from the equation matrix to the solver_matrix_idx'th solver matrix
  END TYPE EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE

  !>Contains information on the mapping from the equations rows in an equations set to the solver rows 
  TYPE EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_SOLVER_ROWS !<The number of solver rows this equations row is mapped to
    INTEGER(INTG), ALLOCATABLE :: SOLVER_ROWS(:) !<SOLVER_ROWS(i). The solver row number for the i'th solver row that this row is mapped to
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The coupling coefficients for the i'th solver row that this row is mapped to.
  END TYPE EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE
  
  !>Contains information on the mappings from an equations set to the solver matrices
  TYPE EQUATIONS_SET_TO_SOLVER_MAP_TYPE
    INTEGER(INTG) :: EQUATIONS_SET_INDEX !<The index of the equations set for these mappings
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mappings
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations in this equations set
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_CONDITIONS !<The number of interface conditions affecting this equations set.
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE_TYPE), ALLOCATABLE :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(:) !<EQUATIONS_TO_SOLVER_MATRIX_MAPS_INTERFACE(interface_condition_idx). Information on the interface_condition_idx'th interface condition affecting this equations set
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM_TYPE), ALLOCATABLE :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(:) !<EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx). The mappings from the equations matrices in this equation set to the solver_matrix_idx'th solver_matrix
    TYPE(EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM_TYPE), ALLOCATABLE :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(:) !<EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(equations_matrix_idx). The mappings from the equation_matrix_idx'th equations matrix in this equation set to the solver_matrices.
    TYPE(JACOBIAN_TO_SOLVER_MAP_PTR_TYPE), ALLOCATABLE :: EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM(:) !<EQUATIONS_TO_SOLVER_MATRIX_MAPS_JM(jacobian_idx). The mappings from the jacobian_idx'th Jacobian matrix in this equation set to the solver_matrices.
    TYPE(EQUATIONS_ROW_TO_SOLVER_ROWS_MAP_TYPE), ALLOCATABLE :: EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(:) !<EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_idx). The mappings from the equations_row_idx'th equations row to the solver matrices rows.
  END TYPE EQUATIONS_SET_TO_SOLVER_MAP_TYPE

  !>Contains information on the interface to solver matrix mappings when indexing by solver matrix number
  TYPE INTERFACE_TO_SOLVER_MATRIX_MAPS_SM_TYPE
    INTEGER(INTG) :: SOLVER_MATRIX_NUMBER !<The number of the solver matrix for these mappings

    INTEGER(INTG) :: LAGRANGE_VARIABLE_TYPE !<LThe variable type for the Lagrange variable that is mapped to the solver matrix.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: LAGRANGE_VARIABLE !<A pointer to the Lagrange variable that is mapped to the solver matrix.
    TYPE(VARIABLE_TO_SOLVER_COL_MAP_TYPE) :: LAGRANGE_VARIABLE_TO_SOLVER_COL_MAP !<The mappings from the Lagrange variable dofs to the solver dofs.
    INTEGER(INTG) :: NUMBER_OF_DEPENDENT_VARIABLES !<The number of dependent variables mapped to this solver matrix
    INTEGER(INTG), ALLOCATABLE :: DEPENDENT_VARIABLE_TYPES(:) !<DEPENDENT_VARIABLE_TYPES(variable_idx). The variable type for the variable_idx'th dependent variable that is mapped to the solver matrix.
    TYPE(FIELD_VARIABLE_PTR_TYPE), ALLOCATABLE :: DEPENDENT_VARIABLES(:) !<DEPENDENT_VARIABLES(variable_idx). DEPENDENT_VARIABLES(variable_idx)%PTR points to the variable_idx'th dependent variable that is mapped to the solver matrix.
    TYPE(VARIABLE_TO_SOLVER_COL_MAP_TYPE), ALLOCATABLE :: DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(:) !<DEPENDENT_VARIABLE_TO_SOLVER_COL_MAPS(interface_matrix_idx). The mappings from the dependent variable dofs to the solver dofs for the interface_matrix_idx'th interface matrix dependent variable in the interface condition that is mapped to the solver matrix.
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_MATRICES !<The number of interface matrices mapped to this solver matrix
    TYPE(INTERFACE_TO_SOLVER_MAPS_PTR_TYPE), ALLOCATABLE :: INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(:) !<INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS(interface_matrix_idx). The maps from the interface)matrix_idx'th interface matrix to solver matrix
    TYPE(EQUATIONS_COL_TO_SOLVER_COLS_MAP_TYPE), ALLOCATABLE :: INTERFACE_COL_TO_SOLVER_COLS_MAP(:) !<EQUATIONS_COL_SOLVER_COL_MAP(column_idx). The mapping from the column_idx'th column of this interface matrix to the solver matrix columns.
  END TYPE INTERFACE_TO_SOLVER_MATRIX_MAPS_SM_TYPE

  !>Contains information on the mapping from an interface condition column to a solver row.
  TYPE INTERFACE_ROW_TO_SOLVER_ROWS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_SOLVER_ROWS !<The number of solver rows that the interface condition row in mapped to (can be either 0 or 1 at he moment.)
    INTEGER(INTG) :: SOLVER_ROW !<The solver row that the interface condition row is mapped to
    REAL(DP) :: COUPLING_COEFFICIENT !<The coupling coefficient for the mapping.
  END TYPE INTERFACE_ROW_TO_SOLVER_ROWS_MAP_TYPE
  
  !>Contains information on the interface to solver matrix mappings when indexing by interface matrix number
  TYPE INTERFACE_TO_SOLVER_MATRIX_MAPS_IM_TYPE
    INTEGER(INTG) :: INTERFACE_MATRIX_NUMBER !<The number of the interface  matrix for these mappings
    INTEGER(INTG) :: NUMBER_OF_SOLVER_MATRICES !<The number of solver matrices mapped to this interface matrix
    TYPE(INTERFACE_TO_SOLVER_MAPS_PTR_TYPE), ALLOCATABLE :: INTERFACE_TO_SOLVER_MATRIX_MAPS(:) !<EQUATIONS_TO_SOLVER_MATRIX_MAPS(solver_matrix_idx). The maps from the equation matrix to the solver_matrix_idx'th solver matrix
    TYPE(INTERFACE_ROW_TO_SOLVER_ROWS_MAP_TYPE), ALLOCATABLE :: INTERFACE_ROW_TO_SOLVER_ROWS_MAP(:) !<INTERFACE_ROW_TO_SOLVER_ROW_MAP(interface_row_idx). The mapping from the interface_row_idx'th interface matrix to a solver row.
  END TYPE INTERFACE_TO_SOLVER_MATRIX_MAPS_IM_TYPE

  TYPE INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS_TYPE
    INTEGER(INTG) :: EQUATIONS_SET_INDEX
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG) :: INTERFACE_MATRIX_INDEX
  END TYPE INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS_TYPE
  
  !>Contains information on the mapping from an interface condition column to a solver row.
  TYPE INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_SOLVER_ROWS !<The number of solver rows that the interface condition column in mapped to (can be either 0 or 1 at he moment.)
    INTEGER(INTG) :: SOLVER_ROW !<The solver row that the interface condition column is mapped to
    REAL(DP) :: COUPLING_COEFFICIENT !<The coupling coefficient for the mapping.
  END TYPE INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP_TYPE
  
  !>Contains information on the mappings from an interface condition to the solver matrices
  TYPE INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE
    INTEGER(INTG) :: INTERFACE_CONDITION_INDEX !<The index of the interface condition for these mappings
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mappings
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations in this interface condition
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS_SETS !<The number of equations sets that the interface condition affects
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS_TYPE), ALLOCATABLE :: INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(:) !<INTERFACE_TO_SOLVER_MATRIX_MAPS_EQUATIONS(equations_set_idx). The equations set information of the equations_set_idx'th equations set that the interface condition affects.
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_SM_TYPE), ALLOCATABLE :: INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(:) !<INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx). The mappings from the interface matrices in this interface condition to the solver_matrix_idx'th solver_matrix
    TYPE(INTERFACE_TO_SOLVER_MATRIX_MAPS_IM_TYPE), ALLOCATABLE :: INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(:) !<INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx). The mappings from the interface_matrix_idx'th interface matrix in this interface condition to the solver_matrices.
    TYPE(INTERFACE_COLUMN_TO_SOLVER_ROWS_MAP_TYPE), ALLOCATABLE :: INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(:) !<INTERFACE_COLUMN_TO_SOLVER_ROW_MAP(interface_column_idx). The mapping from the interface_column_idx'th interface column to a solver row.
  END TYPE INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE
  
 !>Contains information about the mapping from a solver matrix column to dynamic equations matrices and variables
  TYPE SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES !<The number of dynamic equations matrices the solver column is mapped to in this equations set
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_MATRIX_NUMBERS(:) !<EQUATIONS_MATRIX_NUMBERS(i). The i'th equations matrix number in theequations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_COL_NUMBERS(:) !<EQUATIONS_COL_NUMBERS(i). The i'th equations column number in the equation set the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The i'th coupling coefficient for solver column mapping
  END TYPE SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP_TYPE

  !>Contains information about the mapping from a solver matrix column to static equations matrices and variables
  TYPE SOLVER_COL_TO_STATIC_EQUATIONS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_LINEAR_EQUATIONS_MATRICES !<The number of linear equations matrices the solver column is mapped to in this equations set
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_MATRIX_NUMBERS(:) !<EQUATIONS_MATRIX_NUMBERS(i). The i'th equations matrix number in theequations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_COL_NUMBERS(:) !<EQUATIONS_COL_NUMBERS(i). The i'th equations column number in the equation set the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The i'th coupling coefficient for solver column mapping
!!TODO: Maybe split this into a linear and a nonlinear part? The only problem is that would probably use about the same memory???
    INTEGER(INTG) :: JACOBIAN_COL_NUMBER !<The Jacbian column number in the equations set that the solver column is mapped to.
    REAL(DP) :: JACOBIAN_COUPLING_COEFFICIENT !<The coupling coefficient for the solver column to Jacobian column mapping.
  END TYPE SOLVER_COL_TO_STATIC_EQUATIONS_MAP_TYPE

  !>Contains information about mapping the solver dof to the field variable dofs in the equations set.
  TYPE SOLVER_DOF_TO_VARIABLE_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS !<The number of equations this solver dof is mapped to
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_TYPES(:) !<EQUATION_TYPES(equation_idx). The type of equation of the equation_idx'th dof that the solver dof is mapped to.
    INTEGER(INTG), ALLOCATABLE :: EQUATIONS_INDICES(:) !<EQUATIONS_INDICES(equation_idx). The index of either the equations set or interface condition of the equation_idx'th dof that this solver dof is mapped to.
    TYPE(FIELD_VARIABLE_PTR_TYPE), ALLOCATABLE :: VARIABLE(:) !<VARIABLE(equation_idx)%PTR is a pointer to the field variable that the solver dof is mapped to in the equation_idx'th equation
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_DOF(:) !<VARIABLE_DOF(equation_idx). The variable dof number that the solver dof is mapped to in the equation_idx'th equation
    REAL(DP), ALLOCATABLE :: VARIABLE_COEFFICIENT(:) !<VARIABLE_COEFFICIENT(equation_idx). The mulitplicative coefficient for the mapping to the equation_idx'th equations
    REAL(DP), ALLOCATABLE :: ADDITIVE_CONSTANT(:) !<ADDITIVE_CONSTANT(equation_idx). The additive constant for the mapping to the equation_idx'th equations
  END TYPE SOLVER_DOF_TO_VARIABLE_MAP_TYPE
  
  !>Contains information about the mappings from a solver matrix to the equations in an equations set
  TYPE SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations in the equations set that these columns map to.
    LOGICAL :: HAVE_DYNAMIC !<Is .TRUE. if there are any dynamic equations in the map.
    LOGICAL :: HAVE_STATIC !<Is .TRUE. if there are any static equations in the map.
    TYPE(SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAP_TYPE), ALLOCATABLE :: SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(:) !<SOLVER_COL_TO_DYNAMIC_EQUATIONS_MAPS(col_idx). The mappings from the col_idx'th column of the solver matrix to the dynamic equations in the equations set.
    TYPE(SOLVER_COL_TO_STATIC_EQUATIONS_MAP_TYPE), ALLOCATABLE :: SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(:) !<SOLVER_COL_TO_STATIC_EQUATIONS_MAPS(col_idx). The mappings from the col_idx'th column of the solver matrix to the static equations in the equations set.
  END TYPE SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE
  
  !>Contains information about the mapping from a solver matrix column to interface equations matrices and variables
  TYPE SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_MATRICES !<The number of interface matrices the solver column is mapped to in this interface condition
    INTEGER(INTG), ALLOCATABLE :: INTERFACE_MATRIX_NUMBERS(:) !<INTERFACE_MATRIX_NUMBERS(i). The i'th interface matrix number in the interface equations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: INTERFACE_COL_NUMBERS(:) !<INTERFACE_COL_NUMBERS(i). The i'th interface interface column number in the interface condition the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The i'th coupling coefficient for solver column mapping
  END TYPE SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP_TYPE

  !>Contains information about the mappings from a solver matrix to the equations in an equations set
  TYPE SOLVER_COL_TO_INTERFACE_MAP_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations in the interface conditionthat these columns map to.
    TYPE(SOLVER_COL_TO_INTERFACE_EQUATIONS_MAP_TYPE), ALLOCATABLE :: SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS(:) !<SOLVER_COL_TO_INTERFACE_EQUATIONS_MAPS(col_idx). The mappings from the col_idx'th column of the solver matrix to the interface equations in the interface condition.
  END TYPE SOLVER_COL_TO_INTERFACE_MAP_TYPE
  
  !>Contains information on the mappings from a solver matrix to equations sets
  TYPE SOLVER_COL_TO_EQUATIONS_MAPS_TYPE
    INTEGER(INTG) :: SOLVER_MATRIX_NUMBER !<The number of this solver matrix
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping for this solver matrix mapping
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix being mappind
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in this distributed solver matrix
    TYPE(SOLVER_COL_TO_EQUATIONS_SET_MAP_TYPE), ALLOCATABLE :: SOLVER_COL_TO_EQUATIONS_SET_MAPS(:) !<SOLVER_TO_EQUATIONS_SET_MAP(equations_set_idx). The solver columns to equations matrix maps for the equations_set_idx'th equations set.
    TYPE(SOLVER_COL_TO_INTERFACE_MAP_TYPE), ALLOCATABLE :: SOLVER_COL_TO_INTERFACE_MAPS(:) !<SOLVER_COL_TO_INTERFACE_MAPS(interface_condition_idx). The solver columns to interface matrix maps for the interface_condition_idx'th interface condition.
    INTEGER(INTG) :: NUMBER_OF_DOFS !<The number of local dofs (excluding ghost values) in the solver vector associated with this solver matrix
    INTEGER(INTG) :: TOTAL_NUMBER_OF_DOFS !<The number of local dofs (including ghost values) in the solver vector associated with this solver matrix
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_DOFS !<The number of global dofs in the solver vector associated with this solver matrix.
!!TODO: should this be index by solver dof rather than column???
    TYPE(SOLVER_DOF_TO_VARIABLE_MAP_TYPE), ALLOCATABLE :: SOLVER_DOF_TO_VARIABLE_MAPS(:) !<SOLVER_DOF_TO_EQUATIONS_MAPS(dof_idx). The mappings from the dof_idx'th solver dof to the field variables in the equations set.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOFS_MAPPING !<The domain mapping for solver matrix column dofs
  END TYPE SOLVER_COL_TO_EQUATIONS_MAPS_TYPE

  !>Contains information on the mappings from a solver row to the equations.
  TYPE SOLVER_ROW_TO_EQUATIONS_MAPS_TYPE
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS_SETS !<The number of equations sets the solver row is mapped to. If the rows are interface rows then this will be zero.
    INTEGER(INTG) :: INTERFACE_CONDITION_INDEX !<The interface condition index that the solver row is mapped to. If the rows are from an equations set then this will be zero.
     INTEGER(INTG), ALLOCATABLE :: EQUATIONS_INDEX(:) !<EQUATIONS_INDEX(i). If the rows are equations set rows then this is the index of the equations set that the i'th equations set row that the solver row is mapped to. If the rows are interface rows this will not be allocated.
    INTEGER(INTG), ALLOCATABLE :: ROWCOL_NUMBER(:) !<ROWCOL_NUMBER(i). If the row are equations set rows the i'th row number that the solver row is mapped to. If the rows are interface rows then the i'th column number that the solver row is mapped to.
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). The i'th coupling coefficient for the solver row to equations mapping
  END TYPE SOLVER_ROW_TO_EQUATIONS_MAPS_TYPE
  
  !>Contains information about the cached create values for a solver mapping
  TYPE SOLVER_MAPPING_CREATE_VALUES_CACHE_TYPE
    TYPE(LIST_PTR_TYPE), POINTER :: EQUATIONS_VARIABLE_LIST(:) !<EQUATIONS_VARIABLES_LIST(solver_matrix_idx). The list of equations set variables in the solver mapping.
    INTEGER, ALLOCATABLE :: DYNAMIC_VARIABLE_TYPE(:) !<DYNAMIC_VARIABLE_TYPE(equations_set_idx). The variable type that is mapped to the dynamic matrices for the equations_set_idx'th equations set.
    INTEGER(INTG), ALLOCATABLE :: MATRIX_VARIABLE_TYPES(:,:,:) !<MATRIX_VARIABLE_TYPES(0:..,equations_set_idx,matrix_idx). The list of matrix variable types in the equations_set_idx'th equations set for the matrix_idx'th solver matrix. MATRIX_VARIABLE_TYPES(0,equations_set_idx,matrix_idx) is the number of variable types in the equations_set_idx'th equations set mapped to the matrix_idx'th solver matrix and MATRIX_VARIABLE_TYPES(1..,equations_set_idx,matrix_idx) is the list of the variable types in the equations set.
    INTEGER(INTG), ALLOCATABLE :: RESIDUAL_VARIABLE_TYPES(:,:) !<RESIDUAL_VARIABLE_TYPES(0:..,equations_set_idx). The list of residual variable types in the equations_set_idx'th equations set. RESIDUAL_VARIABLE_TYPES(0,equations_set_idx) is the number of variable types in the equations_set_idx'th equations set and RESIDUAL_VARIABLE_TYPES(1..,equations_set_idx) is the list of the variable types in the equations set.
    INTEGER(INTG), ALLOCATABLE :: RHS_VARIABLE_TYPE(:) !<RHS_VARIABLE_TYPE(equations_set_idx). The variable type that is mapped to the solution RHS for the equations_set_idx'th equations set
    INTEGER, ALLOCATABLE :: SOURCE_VARIABLE_TYPE(:) !<SOURCE_VARIABLE_TYPE(equations_set_idx). The source variable type that is mapped to the RHS for the equations_set_idx'th equations set.
    TYPE(LIST_PTR_TYPE), POINTER :: INTERFACE_VARIABLE_LIST(:) !<INTERFACE_VARIABLES_LIST(solver_matrix_idx). The list of interface condition variables in the solver mapping for the solver_matrix idx'th solver matrix.
    TYPE(LIST_PTR_TYPE), POINTER :: INTERFACE_INDICES(:) !<INTERFACE_INDICES(equations_set_idx). The list of interface condition indices in the equations_set_idx'th equations set.
  END TYPE SOLVER_MAPPING_CREATE_VALUES_CACHE_TYPE

  TYPE SOLVER_MAPPING_VARIABLE_TYPE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    INTEGER(INTG) :: VARIABLE_TYPE
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS
    INTEGER(INTG), ALLOCATABLE :: EQUATION_TYPES(:)
    INTEGER(INTG), ALLOCATABLE :: EQUATION_INDICES(:)
  END TYPE SOLVER_MAPPING_VARIABLE_TYPE

  !>Contains information on the variables involved in a solver matrix
  TYPE SOLVER_MAPPING_VARIABLES_TYPE
    INTEGER(INTG) :: NUMBER_OF_VARIABLES !<The number of variables involved in the solver matrix.
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE), ALLOCATABLE :: VARIABLES(:) !<VARIABLES(variable_idx). The variable information for the variable_idx'th variable involved in the solver matrix.
  END TYPE SOLVER_MAPPING_VARIABLES_TYPE
  
  !>Contains information on the solver mapping between the global equation sets and the solver matrices.
  TYPE SOLVER_MAPPING_TYPE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations for this mapping.
    LOGICAL :: SOLVER_MAPPING_FINISHED !<Is .TRUE. if the solution mapping has finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_SOLVER_MATRICES !<The number of solution matrices in this mapping.
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of (local) rows in the solver matrices
    INTEGER(INTG) :: NUMBER_OF_GLOBAL_ROWS !<The number of global rows in the solver matrices
    INTEGER(INTG) :: NUMBER_OF_EQUATIONS_SETS !<The number of equations sets in the solution mapping.
    TYPE(EQUATIONS_SET_PTR_TYPE), ALLOCATABLE :: EQUATIONS_SETS(:) !<The list of equations sets that are in this solution mapping
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE), ALLOCATABLE :: EQUATIONS_SET_TO_SOLVER_MAP(:) !<EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx). The mapping from the equations_set_idx'th equations set to the solver matrices.
    INTEGER(INTG) :: NUMBER_OF_INTERFACE_CONDITIONS !<The number of interface conditions in the solution mapping.
    TYPE(INTERFACE_CONDITION_PTR_TYPE), ALLOCATABLE :: INTERFACE_CONDITIONS(:) !<The list of interface conditions that are in this
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE), ALLOCATABLE :: INTERFACE_CONDITION_TO_SOLVER_MAP(:) !<INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx). The mapping from the interface_condition_idx'th interface condition to the solver matrices.
    TYPE(SOLVER_MAPPING_VARIABLES_TYPE), ALLOCATABLE :: VARIABLES_LIST(:) !<VARIABLES_LIST(solver_matrix_idx). The list of variable for the solver_matrix_idx'th solver matrix.
    TYPE(SOLVER_COL_TO_EQUATIONS_MAPS_TYPE), ALLOCATABLE :: SOLVER_COL_TO_EQUATIONS_COLS_MAP(:) !<SOLVER_TO_EQUATIONS_SET_MAP(solver_matrix_idx). The mapping from the solver_matrix_idx'th solver matrix to the equations set. 
    TYPE(SOLVER_ROW_TO_EQUATIONS_MAPS_TYPE), ALLOCATABLE :: SOLVER_ROW_TO_EQUATIONS_ROWS_MAP(:) !<SOLVER_ROW_TO_EQUATIONS_SET_MAPS(local_row_idx). The mappings from the local_row_idx'th solver row to the equations set rows.
    !LOGICAL :: HAVE_JACOBIAN !<Is .TRUE. if the Jacobian exists for nonlinear problems.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOFS_MAPPING !<The domain mapping for the solver rows.
    TYPE(SOlVER_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<The create values cache for the solver mapping
  END TYPE SOLVER_MAPPING_TYPE

  !
  !================================================================================================================================
  !
  ! History types
  !

  !>Contains information about a history file for a control loop. \see OPENCMISS::CMISSHistoryType
  TYPE HISTORY_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop for the history file
    LOGICAL :: HISTORY_FINISHED !<Is .TRUE. if the history file has finished being created, .FALSE. if not.
    INTEGER(INTG) :: FILE_FORMAT !<The format of the history file \see HISTORY_ROUTINES_FileFormatTypes,HISTORY_ROUTINES
    TYPE(VARYING_STRING) :: FILENAME !<The file name of the history file
    INTEGER(INTG) :: UNIT_NUMBER !<The unit number of the history file.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field that will be read/written in the history file
  END TYPE HISTORY_TYPE
  
  !
  !================================================================================================================================
  !
  ! Control types
  !
  
  !>Contains information on a simple (execute once) control loop
  TYPE CONTROL_LOOP_SIMPLE_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
  END TYPE CONTROL_LOOP_SIMPLE_TYPE

  !>Contains information on a fixed iteration control loop
  TYPE CONTROL_LOOP_FIXED_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    INTEGER(INTG) :: ITERATION_NUMBER
    INTEGER(INTG) :: START_ITERATION
    INTEGER(INTG) :: STOP_ITERATION
    INTEGER(INTG) :: ITERATION_INCREMENT
  END TYPE CONTROL_LOOP_FIXED_TYPE

  !>Contains information on a time iteration control loop
  TYPE CONTROL_LOOP_TIME_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    INTEGER(INTG) :: ITERATION_NUMBER
    INTEGER(INTG) :: GLOBAL_ITERATION_NUMBER
! sebk: is thei usefull?
    INTEGER(INTG) :: OUTPUT_NUMBER
    INTEGER(INTG) :: INPUT_NUMBER
    REAL(DP) :: CURRENT_TIME
    REAL(DP) :: START_TIME
    REAL(DP) :: STOP_TIME
    REAL(DP) :: TIME_INCREMENT
  END TYPE CONTROL_LOOP_TIME_TYPE

  !>Contains information on a do-while control loop
  TYPE CONTROL_LOOP_WHILE_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    INTEGER(INTG) :: ITERATION_NUMBER
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ITERATIONS
    LOGICAL :: CONTINUE_LOOP
  END TYPE CONTROL_LOOP_WHILE_TYPE

  !>Contains information on a load-increment control loop
  TYPE CONTROL_LOOP_LOAD_INCREMENT_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    INTEGER(INTG) :: ITERATION_NUMBER
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ITERATIONS
  END TYPE CONTROL_LOOP_LOAD_INCREMENT_TYPE

  !>A buffer type to allow for an array of pointers to a CONTROL_LOOP_TYPE \see TYPES::CONTROL_LOOP_TYPE
  TYPE CONTROL_LOOP_PTR_TYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: PTR !<The pointer to the control loop
  END TYPE CONTROL_LOOP_PTR_TYPE

  !>Contains information on a control loop. \see OPENCMISS::CMISSControlLoopType
  TYPE CONTROL_LOOP_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer back to the problem for the control loop
    TYPE(CONTROL_LOOP_TYPE), POINTER :: PARENT_LOOP !<A pointer back to the parent control loop if this is a sub loop
    LOGICAL :: CONTROL_LOOP_FINISHED !<Is .TRUE. if the problem control has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: LABEL !<A user defined label for the control loop.
    
    INTEGER(INTG) :: LOOP_TYPE !<The type of control loop \see PROBLEM_CONSTANTS_ControlTypes,PROBLEM_CONSTANTS
    INTEGER(INTG) :: CONTROL_LOOP_LEVEL !<The level of the control loop
    INTEGER(INTG) :: SUB_LOOP_INDEX !<The position of this loop in the sub loops of the parent loop if this is a sub loop.
    INTEGER(INTG) :: OUTPUT_TYPE !<The output type of the control loop \see CONTROL_LOOP_ROUTINES_OutputTypes,CONTROL_LOOP_ROUTINES
    
    TYPE(CONTROL_LOOP_SIMPLE_TYPE), POINTER :: SIMPLE_LOOP !<A pointer to the simple loop information
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: FIXED_LOOP !<A pointer to the fixed loop information
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP !<A pointer to the time loop information
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: WHILE_LOOP !<A pointer to the while loop information
    TYPE(CONTROL_LOOP_LOAD_INCREMENT_TYPE), POINTER :: LOAD_INCREMENT_LOOP !<A pointer to the load increment loop information
    
    INTEGER(INTG) :: NUMBER_OF_SUB_LOOPS !<The number of control loops below this loop
    TYPE(CONTROL_LOOP_PTR_TYPE), ALLOCATABLE :: SUB_LOOPS(:) !<A array of pointers to the loops below this loop.
    
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers for this control loop
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to the history file for this control loop.
  END TYPE CONTROL_LOOP_TYPE  
  
  !
  !================================================================================================================================
  !
  ! Problem types
  !

  TYPE PROBLEM_SETUP_TYPE
    INTEGER(INTG) :: SETUP_TYPE !<The setup type \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
    INTEGER(INTG) :: ACTION_TYPE !<The action type \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
  END TYPE PROBLEM_SETUP_TYPE
  
  !>Contains information for a problem. \see OPENCMISS::CMISSProblemType
  TYPE PROBLEM_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the problem. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the problem in the list of problems.
    LOGICAL :: PROBLEM_FINISHED !<Is .TRUE. if the problem has finished being created, .FALSE. if not.
    TYPE(PROBLEMS_TYPE), POINTER :: PROBLEMS !<A pointer to the problems for this problem.
    
    INTEGER(INTG) :: CLASS !<The problem specification class identifier
    INTEGER(INTG) :: TYPE !<The problem specification type identifier
    INTEGER(INTG) :: SUBTYPE !<The problem specification subtype identifier
    
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop information for the problem.
  END TYPE PROBLEM_TYPE
  
  !>A buffer type to allow for an array of pointers to a PROBLEM_TYPE \see TYPES:PROBLEM_TYPE
  TYPE PROBLEM_PTR_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PTR !<The pointer to the problem.
  END TYPE PROBLEM_PTR_TYPE
       
  !>Contains information on the problems defined.
  TYPE PROBLEMS_TYPE
    INTEGER(INTG) :: NUMBER_OF_PROBLEMS !<The number of problems defined.
    TYPE(PROBLEM_PTR_TYPE), POINTER :: PROBLEMS(:) !<The array of pointers to the problems.
  END TYPE PROBLEMS_TYPE

  !
  !================================================================================================================================
  !
  ! Region types
  
  !>A buffer type to allow for an array of pointers to a REGION_TYPE.
  TYPE REGION_PTR_TYPE
    TYPE(REGION_TYPE), POINTER :: PTR !<The pointer to the region.
  END TYPE REGION_PTR_TYPE
     
  !>Contains information for a region. \see OPENCMISS::CMISSRegionType
  TYPE REGION_TYPE 
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the region. The user number must be unique.
    LOGICAL :: REGION_FINISHED !<Is .TRUE. if the region has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: LABEL !<A user defined label for the region.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system used by the region.
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS  !<A pointer to the data points defined on the region.          
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes defined on the region.
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes defined on the region.
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<A pointer to the generated meshes defined on the region.
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointer to the fields defined on the region.
    TYPE(EQUATIONS_SETS_TYPE), POINTER :: EQUATIONS_SETS !<A pointer to the equation sets defined on the region.
    TYPE(CELLML_ENVIRONMENTS_TYPE), POINTER :: CELLML_ENVIRONMENTS !<A pointer to the CellML environments for the region.
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region for the region. If the region has no parent region then it is the global (world) region and PARENT_REGION is NULL.
    INTEGER(INTG) :: NUMBER_OF_SUB_REGIONS !<The number of sub-regions defined for the region.
    TYPE(REGION_PTR_TYPE), POINTER :: SUB_REGIONS(:) !<An array of pointers to the sub-regions defined on the region. \todo make this allocatable
    TYPE(INTERFACES_TYPE), POINTER :: INTERFACES !<A pointer to the interfaces defined on the region.
  END TYPE REGION_TYPE

  !>Contains information about the regions
  TYPE REGIONS_TYPE
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION !<A pointer to the world region
  END TYPE REGIONS_TYPE

  !
  !================================================================================================================================
  !


END MODULE TYPES

