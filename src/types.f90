!> \file
!> $Id: types.f90 28 2007-07-27 08:35:14Z cpb $
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

  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING

  IMPLICIT NONE
  
  !>Contains information for a particular quadrature scheme. \todo Also evaluate the product of the basis functions at gauss points for speed???
  TYPE QUADRATURE_SCHEME_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the quadrature scheme in the list of quadrature schemes for a particular quadrature.
    TYPE(QUADRATURE_TYPE), POINTER :: QUADRATURE !<The pointer back to the quadrature for a particular quadrature scheme
    INTEGER(INTG) :: NUMBER_OF_GAUSS !<The number of gauss points for the quadrature scheme.
    REAL(DP), ALLOCATABLE :: GAUSS_POSITIONS(:,:) !<GAUSS_POSITIONS(nic,ng). The positions in the nic'th xi coordinate of Gauss point ng. Old CMISS name XIG(ni,ng,nb).
    REAL(DP), ALLOCATABLE :: GAUSS_WEIGHTS(:) !<GAUSS_WEIGHTS(ng). The weight applied to Gauss point ng. Old CMISS name WG(ng,nb).
    REAL(DP), ALLOCATABLE :: GAUSS_BASIS_FNS(:,:,:) !<GAUSS_BASIS_FNS(ns,nu,ng). The value of the basis functions evaluated at Gauss point ng for the nu'th derivative of the basis function associated with the ns'th element parameter. Old CMISS name PG(ns,nu,ng,nb)
  END TYPE QUADRATURE_SCHEME_TYPE

  !>A buffer type to allow for an array of pointers to a QUADRATURE_SCHEME_TYPE \see TYPES::QUADRATURE_SCHEME_TYPE
  TYPE QUADRATURE_SCHEME_PTR_TYPE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: PTR !<A pointer to the quadrature scheme
  END TYPE QUADRATURE_SCHEME_PTR_TYPE

  !>Contains information on the quadrature to be used for integrating a basis.
  TYPE QUADRATURE_TYPE
    INTEGER(INTG) :: TYPE !<The type of the quadrature \see BASIS_ROUTINES_QuadratureTypes
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The pointer back to the basis
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_GAUSS_XI(:) !<NUMBER_OF_GAUSS_XI(ni). For standard Gauss schemes the number of Gauss points to be used in the ni'th xi direction.
    INTEGER(INTG) :: GAUSS_ORDER !<For simplex Gauss schemes the order of the Quadrature scheme i.e., the order/dimension of the polynomial that can be integrated.
    TYPE(QUADRATURE_SCHEME_PTR_TYPE), ALLOCATABLE :: QUADRATURE_SCHEME_MAP(:) !<QUADRATURE_SCHEME_MAP(scheme_idx). The pointer map to the defined quadrature schemes. The size of array is given by BASIS_ROUTINES::BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES. If the quadrature scheme is not defined for the particular type then the array element is NULL. \see BASIS_ROUTINES_QuadratureSchemes.
    INTEGER(INTG) :: NUMBER_OF_SCHEMES !<The number of quadrature schemes defined for this quadrature
    TYPE(QUADRATURE_SCHEME_PTR_TYPE), POINTER :: SCHEMES(:) !<SCHEMES(scheme_idx). The array of pointers to the quadrature schemes defined for the basis. scheme_idx must be between 1 and QUADRATURE_TYPE::NUMBER_OF_SCHEMES.
  END TYPE QUADRATURE_TYPE

  !> A buffer type to allow for an array of pointers to a BASIS_TYPE.
  TYPE BASIS_PTR_TYPE
    TYPE(BASIS_TYPE), POINTER :: PTR !<The pointer to the basis.
  END TYPE BASIS_PTR_TYPE

  !> Contains all information about a basis .
  TYPE BASIS_TYPE
!!TODO: Add in different sub types for the different types of bases???
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
    LOGICAL :: DEGENERATE !<Is .TRUE. if the basis is a degenerate basis (i.e., has collpased nodes), .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). The collpased state of the ni'th direction. COLLAPSED_XI can be either XI_COLLAPSED, COLLAPSED_AT_XI0, COLLAPSED_AT_XI1 or NOT_COLLAPSED dependending on whether or not the ni'th direction is collapsed, has a perpendicular Xi collapsed at the xi=0 end of the ni'th direction, has a perpendicular xi collapsed at the xi=1 of the ni'th direction or is not collapsed. NOTE: in old cmiss the value IBT(1,ni) = 5 or 6 was set for the ni that was collapsed. The perpendicular line/face was then stored in IBT(3,ni). For this the quadratic1 and quadratic2 type interpolation types are set on the perpendicular xi direction and the ni direction that is collapsed will have COLLAPSED_XI(ni) set to XI_COLLAPSED. BE CAREFUL WITH THIS WHEN TRANSLATING OLD CMISS CODE. Old CMISS name IBT(1,ni) ???? \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG) :: NUMBER_OF_COLLAPSED_XI !<The number of xi directions in the basis that are collapsed.
    LOGICAL, ALLOCATABLE :: NODE_AT_COLLAPSE(:) !<NODE_AT_COLLAPSE(nn). Is .TRUE. if the nn'th node of the basis is at a collapse, .FALSE. if not.
    !Quadrature
    TYPE(QUADRATURE_TYPE) :: QUADRATURE !<The quadrature schemes for the basis.
    INTEGER(INTG) :: NUMBER_OF_PARTIAL_DERIVATIVES !<The number of paratial derivatives for the basis. Old CMISS name NUT(nbf)
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of local nodes in the basis. Old CMISS name NNT(nbf)
!!TODO: 
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_NODES_XI(:) !<NUMBER_OF_NODES_XI(ni). The number of local nodes in the ni'th direction in the basis. Old CMISS name IBT(2,ni,nb). \todo CHANGE TO NUMBER_OF_NODES_XIC i.e., xi coordinate index not xi???
    INTEGER(INTG) :: NUMBER_OF_ELEMENT_PARAMETERS  !<The number of element parameters in the basis. Old CMISS name NST(nbf). 
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_DERIVATIVES !<The maximum number of derivatives at any node in the basis. Old CMISS name NKT(0,nbf)
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_DERIVATIVES(:) !<NUMBER_OF_DERIVATIVES(nn). The number of derivatives at the nn'th node in the basis. Old CMISS name NKT(nn,nbf).
    INTEGER(INTG), ALLOCATABLE :: NODE_POSITION_INDEX(:,:) !<NODE_POSITION_INDEX(nn,nic). The index of the node position for the nn'th local node in the nic'th coordinate. For Lagrange-Hermite tensor product basis functions: The number of coordinates equals the number of xi directions. Thus if NODE_POSITION_INDEX(nn,:)=1,2,2 then local node nn is the first node in the ni(c)=1 direction, the second node in the ni(c)=2 direction and the second node in the ni(c)=3 direction; For simplex basis functions: The number of coordinates equals the number of xi directions plus one. The index specifies the inverse distance away from the corner/end of that area coordinate. Thus if an element has quadratic interpolation the index will range from 3 (closest to the corner/end of the element that the area coordinate has the value 1) to 1 (closest to the corners/end of the element that the area coordinate has the value 0). If M is the order of the element then in NODE_POSITION_INDEX(nn,:)=1,1,M then that node is apex for the third area coordinate. In general the index values will add up to M+number of xi directions+1 (i.e., subtract one from the indicies to get the standard simplex coordinates. Old CMISS name INP(nn,ni,nb). \see TYPES::BASIS_TYPE::NODE_POSITION_INDEX_INV.
    INTEGER(INTG), ALLOCATABLE :: NODE_POSITION_INDEX_INV(:,:,:,:) !<NODE_POSITION_INDEX_INV(nnc1,nnc2,nnc3,nnc4). The inverse of the node position index for the basis. The NODE_POSITION_INDEX_INV gives the local node number for the node that has node position indices of nnc1 in the 1st ni(c) direction, nnc2 in the 2nd ni(c) direction, nnc3 in the 3rd ni(c) direction and nnc4 in the 4th ni(c) direction. NOTE: that if the basis has less than 4 ni(c) direction the position index is 1. Old CMISS name NNB(inp1,inp2,inp3,nbf). \see TYPES::BASIS_TYPE::NODE_POSITION_INDEX.
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_ORDER_INDEX(:,:,:) !<DERIVATIVE_ORDER_INDEX(nk,nn,0:ni,nbf). The index of the derivative order for the nk'th derivative of the nn'th node in the ni'th direction of the basis. The derivative index is NO_PART_DERIV for zeroth order, FIRST_PART_DERIV for the first order and SECOND_PART_DERIV for the second order derivative. Thus a DERIVATIVE_ORDER_INDEX(nk,nn,1..) of {NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV} indicates that the nk'th derivative of the nn'th node of the basis is the first derivative with respect to the s2 direction. Old CMISS name IDO(nk,nn,1:ni,nbf). \see TYPES::BASIS_TYPE::DERIVATIVE_ORDER_INDEX_INV,CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_ORDER_INDEX_INV(:,:,:,:) !<DERIVATIVE_ORDER_INDEX_INV(nu1,nu2,nu3,nn). The inverse of the derivative order index for the nn'th local node of the basis. DERIVATIVE_ORDER_INDEX_INV gives the derivative number for the nu1 partial derivative in the 1st xi direction, the nu2 partial derivative in the 2nd xi direction and the nu3 partial derivative in the 3rd xi direction. NOTE: local node nn does not carry any derivatives of the requested partial derivative type then DERIVATIVE_ORDER_INDEX_INV will return 0. If the basis has less than 3 xi directions then the nu index is 1. \see TYPES::BASIS_TYPE::DERIVATIVE_ORDER_INDEX
    INTEGER(INTG), ALLOCATABLE :: PARTIAL_DERIVATIVE_INDEX(:,:) !<PARTIAL_DERIVATIVE_INDEX(nk,nn). Gives the partial derivative number (nu) of the nk'th derivative of the nn'th local node for the basis. Old CMISS name IDO(nk,nn,0,nbf).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_PARAMETER_INDEX(:,:) !<ELEMENT_PARAMETER_INDEX(nk,nn). Gives the element parameter number (ns) of the nk'th derivative of the nn'th local node for the basis. Old CMISS name NSB(nk,nn,nbf).
    !Line information
    INTEGER(INTG) :: NUMBER_OF_LOCAL_LINES !<The number of local lines in the basis.
    INTEGER(INTG), ALLOCATABLE :: LOCAL_LINE_XI_DIRECTION(:) !<LOCAL_LINE_XI_DIRECTION(nae). The Xi direction of the nae'th local line for the basis.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_NODES_IN_LOCAL_LINE(:) !<NUMBER_OF_NODES_IN_LOCAL_LINE(nae). The the number of nodes in the nae'th local line for the basis. Old CMISS name NNL(0,nae,nb).
    INTEGER(INTG), ALLOCATABLE :: NODE_NUMBERS_IN_LOCAL_LINE(:,:) !<NODE_NUMBERS_IN_LOCAL_LINE(nnl,nae). The local node numbers (nn) for the nnl'th line node in the nae'th local line for the basis. Old CMISS name NNL(1..,nae,nb).
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVE_NUMBERS_IN_LOCAL_LINE(:,:) !<DERIVATIVES_NUMBERS_IN_LOCAL_LINE(nnl,nae). The derivative numbers (nk) for the nnl'th line node in the nae'th local line for the basis.
    !Face information
!!TODO:
    !Sub-basis information
    TYPE(BASIS_PTR_TYPE), POINTER :: LINE_BASES(:) !<LINE_BASES(nae). The pointer to the basis for the nae'th line for the basis.
    TYPE(BASIS_PTR_TYPE), POINTER :: FACE_BASES(:) !<FACE_BASES(naf). The pointer to the basis for the naf'th face for the basis.
    INTEGER(INTG) :: NUMBER_OF_SUB_BASES !<The number of sub-bases (lines, faces) for the basis.
    TYPE(BASIS_PTR_TYPE), POINTER :: SUB_BASES(:) !<SUB_BASES(sbn). The pointer to the sbn'th sub-basis for the basis.
    TYPE(BASIS_TYPE), POINTER :: PARENT_BASIS !<The pointer to the parent basis for the basis. NOTE: that if the basis is not a sub-basis of another basis this pointer will be NULL. 
  END TYPE BASIS_TYPE

  !>Contains information on a coordinate system. \todo Have a list of coordinate systems and have a pointer in the coordinate_system_type back to the regions that use them.
  TYPE COORDINATE_SYSTEM_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the coordinate. The user number must be unique.
    LOGICAL :: COORDINATE_SYSTEM_FINISHED !<Is .TRUE. if the coordinate system has finished being created, .FALSE. if not.
    INTEGER(INTG) :: TYPE !<The type of coordinate system. Old CMISS name ITYP10(nr). \see COORINDATE_ROUTINES_CoordinateSystemTypes
    INTEGER(INTG) :: RADIAL_INTERPOLATION_TYPE !<The type of radial interpolation type for non-rectangular cartesian systems. Old CMISS name JTYP10(nr). \see COORDINATE_ROUTINES_RadialInterpolations
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS !<The number of dimensions for the coordinate system. Old CMISS name NJT.
    REAL(DP) :: FOCUS !<The focus of the coordinate system for a prolate-spheriodal coordinate system.
    REAL(DP) :: ORIGIN(3) !<ORIGIN(nj). The nj'th component of the origin of the coordinate system wrt the global coordinate system. NOTE: maybe this should be wrt to the parent regions coordinate system - this would then go into the REGION type.
    REAL(DP) :: ORIENTATION(3,3) !<ORIENTATION(nj,mj). he orientation matrix for the orientation of the coordinate system wrt to the global coordinate system. NOTE: maybe this should be wrt to the parent regions coordinate system - this would then go into the REGION type.
  END TYPE COORDINATE_SYSTEM_TYPE

  !> Contains the interpolated point coordinate metrics. Old CMISS name GL,GU,RG.
  TYPE FIELD_INTERPOLATED_POINT_METRICS_TYPE
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<A pointer to the interpolated point.
    INTEGER(INTG) :: NUMBER_OF_X_DIMENSIONS !<The number of X dimensions.
    INTEGER(INTG) :: NUMBER_OF_XI_DIMENSIONS !<The number of Xi dimensions
    REAL(DP), ALLOCATABLE :: GL(:,:) !<GL(mi,ni). Covariant metric tensor. Old CMISS name GL.
    REAL(DP), ALLOCATABLE :: GU(:,:) !<GU(mi,ni). Contravariant metric tensor. Old CMISS name GU.
    REAL(DP), ALLOCATABLE :: DX_DXI(:,:) !<DX_DXI(nj,ni). Rate of change of the X coordinate system wrt the x coordinate system.
    REAL(DP), ALLOCATABLE :: DXI_DX(:,:) !<DXI_DX(ni,nj). Rate of change of the Xi coordinate system wrt the x coordinate system. 
    REAL(DP) :: JACOBIAN !<The Jacobian of the Xi to X coordinate system transformation. Old CMISS name RG.
    INTEGER(INTG) :: JACOBIAN_TYPE !<The type of Jacobian. \see COORDINATE_ROUTINES_JacobianType
  END TYPE FIELD_INTERPOLATED_POINT_METRICS_TYPE
  
  !>Contains the interpolated value (and the derivatives wrt xi) of a field at a point. Old CMISS name XG.
  TYPE FIELD_INTERPOLATED_POINT_TYPE
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<A pointer to the interpolation parameters of the field that is to be interpolated.
    INTEGER(INTG) :: MAX_PARTIAL_DERIVATIVE_INDEX !<The maximum number of partial derivatives that have been allocated for the values component.
    INTEGER(INTG) :: PARTIAL_DERIVATIVE_TYPE !<The type of the partial derivatives that have been interpolated. PARTIAL_DERIVATIVE_TYPE can be either NO_PART_DERIV, FIRST_PART_DERIV or SECOND_PART_DERIV depending on wether just the field value, the field value and all first derivatives (including cross derivatives) or the first value and all first and second derivatives have been interpolated.
    REAL(DP), ALLOCATABLE :: VALUES(:,:) !<VALUES(component_idx,nu). The interpolated field components and their partial derivatives.
  END TYPE FIELD_INTERPOLATED_POINT_TYPE
  
  !>Contains the parameters required to interpolate a field variable within an element. Old CMISS name XE
  TYPE FIELD_INTERPOLATION_PARAMETERS_TYPE
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to be interpolated.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field VARIABLE to be interpolated.
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: BASES(:) !<BASES(component_idx). An array to hold a pointer to the basis (if any) used for interpolating the component_idx'th component of the field variable.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_PARAMETERS(:) !<NUMBER_OF_PARAMETERS(component_idx). The number of interpolation parameters used for interpolating the component_idx'th component of the field variable.
    REAL(DP), ALLOCATABLE :: PARAMETERS(:,:) !<PARAMETERS(ns,component_idx). The interpolation parameters used for interpolating the component_idx'th component of the field variable.
  END TYPE FIELD_INTERPOLATION_PARAMETERS_TYPE
  
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
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_DOF(:) !<VARIABLE_DOF(ny). The variable dof number for the field ny.
    INTEGER(INTG) :: NUMBER_OF_CONSTANT_DOFS !<The number of constant degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_ELEMENT_DOFS !<The number of element based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_NODE_DOFS !<The number of node based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: NUMBER_OF_POINT_DOFS !<The number of point based degrees-of-freedom in the field dofs.
    INTEGER(INTG), ALLOCATABLE :: CONSTANT_DOF2PARAM_MAP(:,:) !<CONSTANT_DOF2PARAM_MAP(i=1..2,nyy). The mapping from constant field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. When i=2 the DOF2PARAM_MAP gives the variable number (nc) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_DOF2PARAM_MAP(:,:) !<ELEMENT_DOF2PARAM_MAP(i=1..3,nyy). The mapping from element based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2PARAM_MAP gives the element number (ne) of the field parameter. When i=2 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. When i=3 the DOF2PARAM_MAP gives the variable number (nc) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.
    INTEGER(INTG), ALLOCATABLE :: NODE_DOF2PARAM_MAP(:,:) !<NODE_DOF2PARAM_MAP(i=1..4,nyy). The mapping from node based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2PARAM_MAP gives the derivative number (nk) of the field parameter. When i=2 the DOF2PARAM_MAP gives the node number (np) of the field parameter. When i=3 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. When i=4 the DOF2PARAM_MAP gives the variable number (nc) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.
    INTEGER(INTG), ALLOCATABLE :: POINT_DOF2PARAM_MAP(:,:) !<POINT_DOF2PARAM_MAP(i=1..3,nyy). The mapping from point based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2PARAM_MAP gives the point number (nq) of the field parameter. When i=2 the DOF2PARAM_MAP gives the component number (nh) of the field parameter. When i=3 the DOF2PARAM_MAP gives the variable number (nc) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOF_TYPE component of this type.  
  END TYPE FIELD_DOF_TO_PARAM_MAP_TYPE

  !>A type to hold the mapping from field parameters (nodes, elements, etc) to field dof numbers for a particular field variable component.
  TYPE FIELD_PARAM_TO_DOF_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_CONSTANT_PARAMETERS !<The number of constant field parameters for this field variable component. Note: this is currently always 1 but is
  !###      included for completness and to allow for multiple constants per field variable component in the future.
    INTEGER(INTG) :: NUMBER_OF_ELEMENT_PARAMETERS !<The number of element based field parameters for this field variable component.
    INTEGER(INTG) :: NUMBER_OF_NODE_PARAMETERS !<The number of node based field parameters for this field variable component.
    INTEGER(INTG) :: MAX_NUMBER_OF_DERIVATIVES !<The maximum number of derivatives for the node parameters for this field variable component. It is the size of the first index of the NODE_PARAM2DOF_MAP component of this type.
    INTEGER(INTG) :: NUMBER_OF_POINT_PARAMETERS !<The number of point based field parameters for this field variable component.
    INTEGER(INTG) :: CONSTANT_PARAM2DOF_MAP(0:1) !<CONSTANT_PARAM2DOF_MAP(nc). The field dof (nc=0) or variable dof (nc=1) number of the constant parameter for this field variable component. 
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_PARAM2DOF_MAP(:,:) !<ELEMENT_PARAM2DOF_MAP(ne,nc). The field dof (nc=0) or variable dof (nc=1) number of the ne'th element based parameter for this field variable component. \todo Allow for multiple element parameters per element.
    INTEGER(INTG), ALLOCATABLE :: NODE_PARAM2DOF_MAP(:,:,:) !<NODE_PARAM2DOF_MAP(nk,np,nc). The field dof (nc=0) or variable dof (nc=1) number of the nk'th derivative of the np'th node based parameter for this field variable component. Note: because the first index of this array is set to the maximum number of derivatives per node this array wastes memory if there are nodes with a smaller number of derivatives than the maximum. \todo Don't allocate too much memory if there are different numbers of derivatives for different nodes.
    INTEGER(INTG), ALLOCATABLE :: POINT_PARAM2DOF_MAP(:,:) !<POINT_PARAM2DOF_MAP(nq,nc). The field dof (nc=0) or variable dof (nc=1) number of nq'th point based parameter for this field variable component.
  END TYPE FIELD_PARAM_TO_DOF_MAP_TYPE

  !>Contains information for a component of a field variable.
  TYPE FIELD_VARIABLE_COMPONENT_TYPE
    INTEGER(INTG) :: COMPONENT_NUMBER !<The number of the field variable component.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !< A pointer to the field variable for this component.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field for this field variable component.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region for this field variable component.
    INTEGER(INTG) :: INTERPOLATION_TYPE !<The interpolation type of the field variable component \see FIELD_ROUTINES_InterpolationTypes
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The mesh component of the field decomposition for this field variable component.
    INTEGER(INTG) :: SCALING_INDEX !<The index into the defined field scalings for this field variable component.
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain of the field decomposition for this field variable component.
    INTEGER(INTG) :: MAX_NUMBER_OF_INTERPOLATION_PARAMETERS !<The maximum number of interpolations parameters in an element for a field variable component.
    TYPE(FIELD_PARAM_TO_DOF_MAP_TYPE) :: PARAM_TO_DOF_MAP !< The mapping of the field parameters to the field dofs for this field variable component.
  END TYPE FIELD_VARIABLE_COMPONENT_TYPE

  !>Contains information for a field variable defined on a field.
  TYPE FIELD_VARIABLE_TYPE
    INTEGER(INTG) :: VARIABLE_NUMBER !<The number of the field variable
    INTEGER(INTG) :: VARIABLE_TYPE !<The type of the field variable. \see FIELD_ROUTINES_VariableTypes 
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field for this field variable.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region for this field variable.
    INTEGER(INTG) :: MAX_NUMBER_OF_INTERPOLATION_PARAMETERS !<The maximum number of interpolation parameters in an element for a field variable. 
    INTEGER(INTG) :: NUMBER_OF_DOFS !<Number of degress of freedom for this field variable. Old CMISS name NYNR(0,0,nc,nr,nx).
    INTEGER(INTG) :: TOTAL_NUMBER_OF_DOFS !<Total number (global) of degrees of freedom for this field variable. Old CMISS name NYNR(0,0,nc,nr,nx).
    INTEGER(INTG) :: GLOBAL_DOF_OFFSET !<The offset of the start of the global dofs for this variable in the list of global field dofs.
    INTEGER(INTG), ALLOCATABLE :: DOF_LIST(:) !<DOF_LIST(i). The list of field dofs in this field variable. Old CMISS name NYNR(1..,0,nc,nr,nx).
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<Domain mapping for this variable. Only allocated if the field is a dependent field. May have to reconsider this as we now have a domain mapping for the field as a whole and for each variable of the field. The variable domain mapping to is allow a global matrix/vector mapping to be obtained from the field variable.
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !<The number of components in the field variable.
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), ALLOCATABLE :: COMPONENTS(:) !<COMPONENTS(component_idx). The array of field variable components.
  END TYPE FIELD_VARIABLE_TYPE
  
  !>A buffer type to allow for an array of pointers to a FIELD_VARIABLE_TYPE.
  TYPE FIELD_VARIABLE_PTR_TYPE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: PTR !<The pointer to the field variable. 
  END TYPE FIELD_VARIABLE_PTR_TYPE

  !>A type to temporarily hold (cache) the user modifiable values which are used to create a field. 
  TYPE FIELD_CREATE_VALUES_CACHE_TYPE
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !<The number of components in the field for each field variable. NOTE: in the future this will need a variable index on it but just allow for the same number of components for each field variable for now. 
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). The cache of the variable type for the given variable_idx of the field. \see FIELD_ROUTINES_VariableTypes
    INTEGER(INTG), ALLOCATABLE :: INTERPOLATION_TYPE(:,:) !<INTERPOLATION_TYPES(component_idx,variable_idx). The cache of the interpolation type for the given component and variable of the field. \see FIELD_ROUTINES_InterpolationTypes
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENT_NUMBER(:,:) !<MESH_COMPONENT_NUMBER(component_idx,varaible_idx). The cache of the mesh component number for the given component and variable of the field.
  END TYPE FIELD_CREATE_VALUES_CACHE_TYPE

  !>The type containing the mappings for the field
  TYPE FIELD_MAPPINGS_TYPE
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE) :: DOF_TO_PARAM_MAP !<The mappings for the field dofs to the field parameters
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The domain mappings for the field dofs i.e., the local to global mpapings etc.
  END TYPE FIELD_MAPPINGS_TYPE

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
    INTEGER(INTG) :: NUMBER_OF_PARAMETER_SETS !<The number of parameter sets that are currently defined on the field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field that these parameter sets are defined on.
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE), POINTER :: SET_TYPE(:) !<SET_TYPE(set_type_idx). A pointer to an array of pointers to the field set types. SET_TYPE(set_type_idx)%PTR is a pointer to the parameter set type for the set_type_idx'th parameter set. set_type_idx can vary from 1 to FIELD_ROUTINES::FIELD_NUMBER_OF_SET_TYPES. The value of the pointer will be NULL if the parameter set corresponding to the set_type_idx'th parameter set has not yet been created for the field.
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE), POINTER :: PARAMETER_SETS(:) !<PARAMETER_SETS(set_type_idx). A pointer to an array of pointers to the parameter sets that have been created on the field. PARAMETER_SET(set_type_idx)%PTR is a pointer to the parameter set type for the set_type_idx'th parameter set that has been created. set_type_idx can vary from 1 to the number of parameter set types that have currently been created for the field i.e., TYPES::FIELD_PARAMETER_SETS_TYPE::NUMBER_OF_PARAMETER_SETS.
  END TYPE FIELD_PARAMETER_SETS_TYPE
  
  !>Contains information for a field defined on a region.
  TYPE FIELD_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the field in the list of fields for a region.
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the field. The user number must be unique.
    LOGICAL :: FIELD_FINISHED !<Is .TRUE. if the field has finished being created, .FALSE. if not.
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointer to the fields for this region.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region for this field.
    INTEGER(INTG) :: TYPE !<The type of the field. \see FIELD_ROUTINES_FieldTypes
    INTEGER(INTG) :: DEPENDENT_TYPE !<The dependent type of the field. \see FIELD_ROUTINES_DependentTypes
    INTEGER(INTG) :: DIMENSION !<The dimension of the field. \see FIELD_ROUTINES_DimensionTypes
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition of the mesh for which the field is defined on.
    INTEGER(INTG) :: NUMBER_OF_VARIABLES !<The number of variable types in the field. Old CMISS name NCT(nr,nx)
    TYPE(FIELD_VARIABLE_PTR_TYPE), ALLOCATABLE :: VARIABLE_TYPE_MAP(:) !<VARIABLE_TYPE_MAP(variable_idx). The map from the available field variable types to the field variable types that are defined for the field. variable_idx varies from 1 to FIELD_ROUTINES::FIELD_NUMBER_OF_VARIABLE_TYPES. If the particular field variable type has not been defined on the field then the VARIABLE_TYPE_MAP will be NULL. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_VARIABLE_TYPE), ALLOCATABLE :: VARIABLES(:) !<VARIABLES(variable_idx) .The array of field variables. 
    TYPE(FIELD_SCALINGS_TYPE) :: SCALINGS !<The scaling parameters for the field
    TYPE(FIELD_MAPPINGS_TYPE) :: MAPPINGS !<The mappings for the field
    TYPE(FIELD_PARAMETER_SETS_TYPE) :: PARAMETER_SETS !<The parameter sets for the field
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field that this field uses. If the field itself is a geometric field then this will be a pointer back to itself.
    TYPE(FIELD_GEOMETRIC_PARAMETERS_TYPE), POINTER :: GEOMETRIC_FIELD_PARAMETERS !<
    TYPE(FIELD_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<The create values cache for the field.
  END TYPE FIELD_TYPE

  !>A buffer type to allow for an array of pointers to a FIELD_TYPE.
  TYPE FIELD_PTR_TYPE
    TYPE(FIELD_TYPE), POINTER :: PTR !<The pointer to the field.  
  END TYPE FIELD_PTR_TYPE

  !>Contains information on the fields defined on a region.
  TYPE FIELDS_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the fields.
    INTEGER(INTG) :: NUMBER_OF_FIELDS !<The number of fields defined on the region.
    TYPE(FIELD_PTR_TYPE), POINTER :: FIELDS(:) !<FIELDS(fields_idx). The array of pointers to the fields.
  END TYPE FIELDS_TYPE

  !>Contains information about a node.
  TYPE NODE_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of node.
    INTEGER(INTG) :: USER_NUMBER !<The user defined number of node.
    TYPE(VARYING_STRING) :: LABEL !<A string label for the node
    REAL(DP), ALLOCATABLE :: INITIAL_POSITION(:) !<INITIAL_POSITION(nj). The initial position of of the node. Used to identify the position of the node before meshing (e.g., Delauny triangulisation) as the geometric field has not been created when the node is created. The actual geometric coordinates for computation using the node should be taken from the geometric field which uses the node.
  END TYPE NODE_TYPE

  !>Contains information on the nodes defined on a region.
  TYPE NODES_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the nodes.
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of nodes defined on the region.
    LOGICAL :: NODES_FINISHED !<Is .TRUE. if the nodes have finished being created, .FALSE. if not.
    TYPE(NODE_TYPE), POINTER :: NODES(:) !<NODES(nodes_idx). A pointer to the nodes. \todo Should this be allocatable?
  END TYPE NODES_TYPE

  !>Contains information on the dofs for a mesh.
  TYPE MESH_DOFS_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh.
    INTEGER(INTG) :: NUMBER_OF_DOFS !<The number of dofs in the mesh.
  END TYPE MESH_DOFS_TYPE

  !>Contains the information for an element in a mesh.
  TYPE MESH_ELEMENT_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global element number in the mesh.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the element.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function for the element.
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:) !<GLOBAL_ELEMENT_NODES(nn). The global node number in the mesh of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: USER_ELEMENT_NODES(:) !<USER_ELEMENT_NODES(nn). The user node number in the mesh of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ADJACENT_ELEMENTS(:) !<NUMBER_OF_ADJACENT_ELEMENTS(-ni:ni). The number of elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent element before the element in the ni'th direction and +ni gives the adjacent element after the element in the ni'th direction. The ni=0 index should be 1 for the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:,:) !<ADJACENT_ELEMENTS(nei,-ni:ni). The local element numbers of the elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent elements before the element in the ni'th direction and +ni gives the adjacent elements after the element in the ni'th direction. The ni=0 index should give the current element number. Old CMISS name NXI(-ni:ni,0:nei,ne)
  END TYPE MESH_ELEMENT_TYPE

  !>Contains the information for the elements of a mesh.
  TYPE MESH_ELEMENTS_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<The pointer to the mesh for the elements information.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !< The number of elements in the mesh.
    LOGICAL :: ELEMENTS_FINISHED !<Is .TRUE. if the mesh elements have finished being created, .FALSE. if not.
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENTS(:) !<ELEMENTS(ne). The pointer to the array of information for the elements of this mesh. ELEMENTS(ne) contains the information for the ne'th global element of the mesh. \todo Should this be allocatable.
  END TYPE MESH_ELEMENTS_TYPE

  !>Contains the topology information for a global node of a mesh.
  TYPE MESH_NODE_TYPE
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global node number in the mesh.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the node.
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES !<The number of global derivatives at the node for the mesh. Old CMISS name NKT(nj,np).
    INTEGER(INTG), ALLOCATABLE :: PARTIAL_DERIVATIVE_INDEX(:) !<PARTIAL_DERIVATIVE_INDEX(nk). The partial derivative index (nu) of the nk'th global derivative for the node. Old CMISS name NUNK(nk,nj,np).
    INTEGER(INTG), ALLOCATABLE :: DOF_INDEX(:) !<DOF_INDEX(nk). The global dof derivative index (ny) in the domain of the nk'th global derivative for the node.
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements surrounding the node in the mesh. Old CMISS name NENP(np,0,0:nr).
    INTEGER(INTG), POINTER :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nep). The global element number of the nep'th element that is surrounding the node. Old CMISS name NENP(np,nep,0:nr). \todo Change this to allocatable.
  END TYPE MESH_NODE_TYPE

  !>Contains the information for the nodes of a mesh.
  TYPE MESH_NODES_TYPE
    TYPE(MESH_TYPE), POINTER :: MESH !<The pointer to the mesh for this nodes information.
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of nodes in the mesh.
    TYPE(MESH_NODE_TYPE), POINTER :: NODES(:) !<NODES(np). The pointer to the array of topology information for the nodes of the mesh. NODES(np) contains the topological information for the np'th global node of the mesh. \todo Should this be allocatable???
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
  
  !>Contains information on the degrees-of-freedom (dofs) for a domain.
  TYPE DOMAIN_DOFS_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain.
    INTEGER(INTG) :: NUMBER_OF_DOFS !<The number of degrees-of-freedom (excluding ghost dofs) in the domain.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_DOFS !<The total number of degrees-of-freedom (including ghost dofs) in the domain.
    INTEGER(INTG), ALLOCATABLE :: DOF_INDEX(:,:) !<DOF_INDEX(i,ny). The index for the ny'th degree-of-freedom. When i=1 DOF_INDEX will give the global derivative number (nk) associated with the dof. When i=2 DOF_INDEX will give the local node number (np) associated with the dof.
  END TYPE DOMAIN_DOFS_TYPE

  !>Contains the information for a line in a domain.
  TYPE DOMAIN_LINE_TYPE
    INTEGER(INTG) :: NUMBER !<The line number in the domain.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function for the line.
    INTEGER(INTG), ALLOCATABLE :: NODES_IN_LINE(:) !<NODES_IN_LINE(nn). The local node number in the domain of the nn'th local node in the line. Old CMISS name NPL(2..5,nj,nl).
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVES_IN_LINE(:,:) !<DERIVATIVES_IN_LINE(nk,nn). The global derivative number of the local derivative nk for the local node nn in the line. Old CMISS name NPL(4..5,nj,nl).
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
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVES_IN_FACE(:,:) !<DERIVATIVES_IN_FACE(nk,nn). The global derivative number of the local derivative nk for the local node nn in the face.
  END TYPE DOMAIN_FACE_TYPE

  !>A buffer type to allow for an array of pointers to a FIELD_VARIABLE_TYPE.
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
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:) !<ELEMENT_NODES(nn). The local node number in the domain of the nn'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_DERIVATIVES(:,:) !<ELEMENT_DERIVATIVES(nk,nn). The global derivative number of the local derivative nk for the local node nn in the element. Old CMISS name NKJE(nk,nn,nj,ne).
  END TYPE DOMAIN_ELEMENT_TYPE
  
  !>Contains the topology information for the elements of a domain.
  TYPE DOMAIN_ELEMENTS_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this elements topology information.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !<The number of elements (excluding ghost elements) in this domain topology.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ELEMENTS !<The total number of elements (including ghost elements) in this domain topology.
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: ELEMENTS(:) !<ELEMENTS(ne). The pointer to the array of topology information for the elements of this domain. ELEMENTS(ne) contains the topological information for the ne'th local elements of the domain. \todo Change this to allocatable???
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS !<The maximum number of element parameters (ns) for all the elements in the domain.
  END TYPE DOMAIN_ELEMENTS_TYPE

  !>Contains the topology information for a local node of a domain.
  TYPE DOMAIN_NODE_TYPE
    INTEGER(INTG) :: LOCAL_NUMBER !<The local node number in the domain.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global node number in the mesh of the local node number in the domain.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the node.
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES !<The number of global derivatives at the node for the domain. Old CMISS name NKT(nj,np)
    INTEGER(INTG), ALLOCATABLE :: PARTIAL_DERIVATIVE_INDEX(:) !<PARTIAL_DERIVATIVE_INDEX(nk). The partial derivative index (nu) of the nk'th global derivative for the node. Old CMISS name NUNK(nk,nj,np).
    INTEGER(INTG), ALLOCATABLE :: DOF_INDEX(:) !<DOF_INDEX(nk). The local dof derivative index (ny) in the domain of the nk'th global derivative for the node.
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements surrounding the node in the domain. Old CMISS name NENP(np,0,0:nr).
    INTEGER(INTG), POINTER :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nep). The local element number of the nep'th element that is surrounding the node. Old CMISS name NENP(np,nep,0:nr). \todo Change this to allocatable.
    INTEGER(INTG) :: NUMBER_OF_NODE_LINES !<The number of lines surrounding the node in the domain.
    INTEGER(INTG), ALLOCATABLE :: NODE_LINES(:) !<NODE_LINES(nlp). The local line number of the nlp'th line that is surrounding the node.
  END TYPE DOMAIN_NODE_TYPE

  !>Contains the topology information for the nodes of a domain
  TYPE DOMAIN_NODES_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this nodes topology information.
    INTEGER(INTG) :: NUMBER_OF_NODES !<The number of nodes (excluding ghost nodes) in this domain topology.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of nodes (including ghost nodes) in this domain topology.
    INTEGER(INTG) :: MAXIMUM_NUMBER_OF_DERIVATIVES !<The maximum number of derivatives over the nodes in this domain topology.
    TYPE(DOMAIN_NODE_TYPE), POINTER :: NODES(:) !<NODES(np). The pointer to the array of topology information for the nodes of this domain. NODES(np) contains the topological information for the np'th local node of the domain. \todo Change this to allocatable???
  END TYPE DOMAIN_NODES_TYPE

  !>Contains the topology information for a domain
  TYPE DOMAIN_TOPOLOGY_TYPE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<The pointer to the domain for this topology information.
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES !<The pointer to the topology information for the nodes of this domain.
    TYPE(DOMAIN_DOFS_TYPE), POINTER :: DOFS !<The pointer to the topology information for the dofs of this domain.
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS !<The pointer to the topology information for the elements of this domain.
    TYPE(DOMAIN_FACES_TYPE), POINTER :: FACES !<The pointer to the topology information for the faces of this domain.
    TYPE(DOMAIN_LINES_TYPE), POINTER :: LINES !<The pointer to the topology information for the lines of this domain.
  END TYPE DOMAIN_TOPOLOGY_TYPE

  !>Contains the information for an adjacent domain for transfering the ghost data of a distributed vector to/from the
  !>current domain.
  TYPE DISTRIBUTED_VECTOR_TRANSFER_TYPE
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR !<The pointer to the distributed vector object for this transfer information.
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

  !>Contains the information for a vector that is distributed across a number of domains.
  TYPE DISTRIBUTED_VECTOR_TYPE
    INTEGER(INTG) :: BASE_TAG_NUMBER !<The base number for the MPI tag numbers that will be used to communicate the distributed vector data amongst the domains. The base tag number can be thought of as the identification number for the distributed vector object.
    LOGICAL :: VECTOR_FINISHED !<!<Is .TRUE. if the distributed vector has finished being created, .FALSE. if not.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The pointer for the domain mapping that identifies how the vector is distributed amongst the domains.
    INTEGER(INTG) :: N !<The size of the distributed vector
    INTEGER(INTG) :: DATA_TYPE !<The type of data for the distributed vector \see DISTRIBUTED_MATRIX_VECTOR_DataTypes 
    INTEGER(INTG) :: DATA_SIZE !<The size of the distributed vector that is held locally by the domain.
    INTEGER(INTG), ALLOCATABLE :: DATA_INTG(:) !<DATA_INTG(i). The integer data for an integer distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    REAL(DP), ALLOCATABLE :: DATA_DP(:) !<DATA_DP(i). The real data for a double precision real distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    REAL(SP), ALLOCATABLE :: DATA_SP(:) !<DATA_SP(i). The real data for a single precision real distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    LOGICAL, ALLOCATABLE :: DATA_L(:) !<DATA_L(i). The logical data for a logical distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain.  
    TYPE(DISTRIBUTED_VECTOR_TRANSFER_TYPE), ALLOCATABLE :: TRANSFERS(:) !<TRANSFERS(adjacent_domain_idx). The transfer information for the adjacent_domain_idx'th adjacent domain to this domain. 
  END TYPE DISTRIBUTED_VECTOR_TYPE

  !>Contains the information for a matrix that is distributed across a number of domains.
  TYPE DISTRIBUTED_MATRIX_TYPE
    INTEGER(INTG) :: BASE_TAG_NUMBER !<The base number for the MPI tag numbers that will be used to communicate the distributed matrix data amongst the domains. The base tag number can be thought of as the identification number for the distributed matrix object.
    LOGICAL :: MATRIX_FINISHED !<Is .TRUE. if the distributed matrix has finished being created, .FALSE. if not.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The pointer for the domain mapping that identifies how the matrix is distributed amongst the domains.
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to store the rows corresponding to this domain.
  END TYPE DISTRIBUTED_MATRIX_TYPE

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
    INTEGER(INTG), ALLOCATABLE :: ROW_INDICES(:) !<ROW_INDICES(i). The row indices for the matrix storage scheme. \see MATRIX_VECTOR_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column indices for the matrix storage scheme. \see MATRIX_VECTOR_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: DATA_INTG(:) !<DATA_INTG(i). The integer data for an integer matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    REAL(SP), ALLOCATABLE :: DATA_SP(:) !<DATA_SP(i). The real data for a single precision matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    REAL(DP), ALLOCATABLE :: DATA_DP(:) !<DATA_DP(i). The real data for a double precision matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    LOGICAL, ALLOCATABLE :: DATA_L(:) !<DATA_L(i). The logical data for a logical matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
  END TYPE MATRIX_TYPE
  
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
    INTEGER(INTG), ALLOCATABLE :: LOCAL_TYPE(:) !<LOCAL_TYPE(domain_idx). The type of local for the domain_idx'th domain for which the global number is mapped to a local number. The types depend on wether the mapped local number in the domain_idx'th domain is an internal, boundary or ghost local number. \see DOMAIN_MAPPINGS_DomainType
  END TYPE DOMAIN_GLOBAL_MAPPING_TYPE

  !>Contains information on the domain mappings (i.e., local and global numberings).
  TYPE DOMAIN_MAPPING_TYPE
    INTEGER(INTG) :: NUMBER_OF_LOCAL !<The number of local numbers in the domain excluding ghost numbers
    INTEGER(INTG) :: TOTAL_NUMBER_OF_LOCAL !<The total number of local numbers in the domain including ghost numbers.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_DOMAIN_LOCAL(:) !<NUMBER_OF_DOMAIN_LOCAL(domain_no). The total number of locals for domain_no'th domain. NOTE: the domain_no goes from 0 to the number of domains-1.
    INTEGER(INTG) :: NUMBER_OF_GLOBAL !<The number of global numbers for this mapping.
    INTEGER(INTG) :: NUMBER_OF_DOMAINS !<The number of domains in this mapping.
    INTEGER(INTG) :: NUMBER_OF_INTERNAL !<The number of internal numbers in this mapping.
    INTEGER(INTG), ALLOCATABLE :: INTERNAL_LIST(:) !<INTERNAL_LIST(i). The list of internal numbers for the mapping. The i'th position gives the i'th local internal number.
    INTEGER(INTG) :: NUMBER_OF_BOUNDARY !<The number of boundary numbers in this mapping.
    INTEGER(INTG), ALLOCATABLE :: BOUNDARY_LIST(:) !<BOUNDARY_LIST(i). The list of boundary numbers for the mapping. The i'th position gives the i'th local boundary number.
    INTEGER(INTG) :: NUMBER_OF_GHOST !<The number of ghost numbers in this mapping.
    INTEGER(INTG), ALLOCATABLE :: GHOST_LIST(:) !<GHOST_LIST(i). The list of ghost numbers for the mapping. The i'th position gives the i'th local ghost number.
    !!Need a dimension index here as domain  variables that map to matrices will have different row and column
    !!mappings in general????
    INTEGER(INTG), ALLOCATABLE :: LOCAL_TO_GLOBAL_MAP(:) !<LOCAL_TO_GLOBAL_MAP(i). The global number for the i'th local number for the mapping.
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), ALLOCATABLE :: GLOBAL_TO_LOCAL_MAP(:) !<GLOBAL_TO_LOCAL_MAP(i). The local information for the i'th global number for the mapping.
    INTEGER(INTG) :: NUMBER_OF_ADJACENT_DOMAINS !<The number of domains that are adjacent to this domain in the mapping.
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_DOMAINS_PTR(:) !<ADJACENT_DOMAINS_PTR(domain_no). The pointer to the list of adjacent domains for domain_no. ADJACENT_DOMAINS_PTR(domain_no) gives the starting position in ADJACENT_DOMAINS_LIST for the first adjacent domain number for domain number domain_no. ADJACENT_DOMAINS_PTR(domain_no+1) gives the last+1 position in ADJACENT_DOMAINS_LIST for the last adjacent domain number for domain number domain_no. NOTE: the index for ADJACENT_DOMAINS_PTR varies from 0 to the number of domains+1.
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
  
  !>Contains the information for a line in a decomposition.
  TYPE DECOMPOSITION_LINE_TYPE
    INTEGER(INTG) :: NUMBER !<The line number in the decomposition.
    INTEGER(INTG) :: XI_DIRECTION !<The Xi direction of the line. Old CMISS name NPL(1,0,nl)
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS !<The number of elements that surround (use) this line.
    INTEGER(INTG), ALLOCATABLE :: SURROUNDING_ELEMENTS(:) !<SURROUNDING_ELEMENTS(nel). The local element number of the nel'th element that surrounds (uses) this line. 
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LINES(:) !<ELEMENT_LINES(nel). The local arc number of the nel'th element that surrounds (uses) this line.
    INTEGER(INTG) :: ADJACENT_LINES(0:1) !<ADJACENT_LINES(0:1). The line number of adjacent lines. ADJACENT_LINES(0) is the line number adjacent in the -xi direction. ADJACENT_LINES(1) is the line number adjacent in the +xi direction. Old CMISS name NPL(2..3,0,nl).
  END TYPE DECOMPOSITION_LINE_TYPE

  !>Contains the topology information for the lines of a decomposition.
  TYPE DECOMPOSITION_LINES_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this lines topology information.
    INTEGER(INTG) :: NUMBER_OF_LINES !<The number of lines in this decomposition topology.
    TYPE(DECOMPOSITION_LINE_TYPE), ALLOCATABLE :: LINES(:) !<LINES(nl). The pointer to the array of topology information for the lines of this decomposition. LINES(nl) contains the topological information for the nl'th local line of the decomposition.
  END TYPE DECOMPOSITION_LINES_TYPE

  !>Contains the information for an element in a decomposition.
  TYPE DECOMPOSITION_ELEMENT_TYPE
    INTEGER(INTG) :: LOCAL_NUMBER !<The local element number in the decomposition.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global element number in the mesh of the local element number in the decomposition.
    INTEGER(INTG) :: USER_NUMBER !<The corresponding user number for the element.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ADJACENT_ELEMENTS(:) !<NUMBER_OF_ADJACENT_ELEMENTS(-ni:ni). The number of elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent element before the element in the ni'th direction and +ni gives the adjacent element after the element in the ni'th direction. The ni=0 index should be 1 for the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:,:) !<ADJACENT_ELEMENTS(nei,-ni:ni). The local element numbers of the elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent elements before the element in the ni'th direction and +ni gives the adjacent elements after the element in the ni'th direction. The ni=0 index should give the current element number. Old CMISS name NXI(-ni:ni,0:nei,ne).
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LINES(:) !<ELEMENT_LINES(nae). The local decomposition line number corresponding to the nae'th local line of the element. Old CMISS name NLL(nae,ne). 
  END TYPE DECOMPOSITION_ELEMENT_TYPE

  !>Contains the topology information for the elements of a decomposition.
  TYPE DECOMPOSITION_ELEMENTS_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this elements topology information.
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ELEMENTS !<The total number of elements in this decomposition topology.
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: ELEMENTS(:) !<ELEMENTS(ne). The pointer to the array of topology information for the elements of this decomposition. ELEMENTS(ne) contains the topological information for the ne'th local element of the decomposition. \todo Change this to allocatable???
  END TYPE DECOMPOSITION_ELEMENTS_TYPE

   !>Contains the topology information for a decomposition
  TYPE DECOMPOSITION_TOPOLOGY_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<The pointer to the decomposition for this topology information.
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS !<The pointer to the topology information for the elements of this decomposition.
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: LINES !<The pointer to the topology information for the lines of this decomposition.
  END TYPE DECOMPOSITION_TOPOLOGY_TYPE

  !>Contains information on the domain decomposition.
  TYPE DECOMPOSITION_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the domain decomposition. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the domain decomposition in the list of domain decompositions for a particular mesh.
    LOGICAL :: DECOMPOSITION_FINISHED !<Is .TRUE. if the decomposition has finished being created, .FALSE. if not.
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS !<A pointer to the decompositions for this decomposition.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh for this decomposition.
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER !<The component number (index) of the mesh component that this decomposition belongs to (i.e., was generated from).
    INTEGER(INTG) :: DECOMPOSITION_TYPE !<The type of the domain decomposition \see MESH_ROUTINES_DecompositionTypes.
    INTEGER(INTG) :: NUMBER_OF_DOMAINS !<The number of domains that this decomposition contains.
    INTEGER(INTG) :: NUMBER_OF_EDGES_CUT !<For automatically calcualted decompositions, the number of edges of the mesh dual graph that were cut for the composition. It provides an indication of the optimally of the automatic decomposition.
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_DOMAIN(:) !<ELEMENT_DOMAIN(ne). The domain number that the ne'th global element is in for the decomposition. Note: the domain numbers start at 0 and go up to the NUMBER_OF_DOMAINS-1.
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the topology for this decomposition.
    TYPE(DOMAIN_PTR_TYPE), POINTER :: DOMAIN(:) !<DOMAIN(mesh_component_idx). A pointer to the domain for mesh component for the domain associated with the computational node. \todo Change this to allocatable???
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

  !>Contains information on a generated regular mesh
  TYPE GENERATED_MESH_REGULAR_TYPE
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    REAL(DP), ALLOCATABLE :: ORIGIN(:) !<ORIGIN(nj). The position of the origin (first) corner of the regular mesh
    REAL(DP), ALLOCATABLE :: MAXIMUM_EXTENT(:) !<MAXIMUM_EXTENT(nj). The extent/size in each nj'th direction of the regular mesh.
    INTEGER(INTG) :: MESH_DIMENSION !<The dimension/number of Xi directions of the regular mesh.
    INTEGER(INTG), ALLOCATABLE :: NUMBER_OF_ELEMENTS_XI(:) !<NUMBER_OF_ELEMENTS_XI(ni). The number of elements in the ni'th Xi direction for the mesh.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The pointer to the basis used in the regular mesh.
  END TYPE GENERATED_MESH_REGULAR_TYPE

  !>Contains information for generated meshes.
  TYPE GENERATED_MESH_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user number of the generated mesh.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the generated mesh.
    INTEGER(INTG) :: TYPE !<The type of generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH !<A pointer to the information for a regular generated mesh. 
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh that is generated
  END TYPE GENERATED_MESH_TYPE
  
  !>Contains information on a mesh defined on a region.
  TYPE MESH_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user number of the mesh. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The corresponding global number for the mesh.
    LOGICAL :: MESH_FINISHED !<Is .TRUE. if the mesh has finished being created, .FALSE. if not.
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes for this mesh.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing this mesh.
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS !<The number of dimensions (Xi directions) for this mesh.
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !<The number of mesh components in this mesh.
    LOGICAL :: MESH_EMBEDDED !<Is .TRUE. if the mesh is embedded in another mesh, .FALSE. if not.
    TYPE(MESH_TYPE), POINTER :: EMBEDDING_MESH !<If this mesh is embedded the pointer to the mesh that this mesh is embedded in. IF the mesh is not embedded the pointer is NULL.
    INTEGER(INTG) :: NUMBER_OF_EMBEDDED_MESHES !<The number of meshes that are embedded in this mesh.
    TYPE(MESH_PTR_TYPE), POINTER :: EMBEDDED_MESHES(:) !<EMBEDDED_MESHES(mesh_idx). A pointer to the mesh_idx'th mesh that is embedded in this mesh.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !<The number of elements in the mesh.
    INTEGER(INTG) :: NUMBER_OF_FACES !<The number of faces in the mesh.
    INTEGER(INTG) :: NUMBER_OF_LINES !<The number of lines in the mesh.
    TYPE(MESH_TOPOLOGY_PTR_TYPE), POINTER :: TOPOLOGY(:) !<TOPOLOGY(mesh_component_idx). A pointer to the topology mesh_component_idx'th mesh component.
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS !<A pointer to the decompositions for this mesh.
  END TYPE MESH_TYPE

  !>A buffer type to allow for an array of pointers to a MESH_TYPE.
  TYPE MESH_PTR_TYPE
    TYPE(MESH_TYPE), POINTER :: PTR !<The pointer to the mesh. 
  END TYPE MESH_PTR_TYPE

  !>Contains information on the meshes defined on a region.
  TYPE MESHES_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region.
    INTEGER(INTG) :: NUMBER_OF_MESHES !<The number of meshes defined on the region.
    TYPE(MESH_PTR_TYPE), POINTER :: MESHES(:) !<MESHES(meshes_idx). The array of pointers to the meshes.
  END TYPE MESHES_TYPE

  !>Contains information on the geometry for a problem
  TYPE PROBLEM_GEOMETRY_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem.
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<The geometric field for this problem.
    TYPE(FIELD_TYPE), POINTER :: FIBRE_FIELD !<The fibre field for this problem if one is defined. If no fibre field is defined the pointer is NULL.
  END TYPE PROBLEM_GEOMETRY_TYPE

  !>Contains information on the materials for the problem.
  TYPE PROBLEM_MATERIALS_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem.
    LOGICAL :: MATERIALS_FINISHED !<Is .TRUE. if the materials for the problem has finished being created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: MATERIAL_FIELD !<A pointer to the material field for the problem if one is defined. If no material field is defined the pointer is NULL.
  END TYPE PROBLEM_MATERIALS_TYPE

  !>Contains information on the fixed conditions for the problem.
  TYPE PROBLEM_FIXED_CONDITIONS_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem.
    LOGICAL :: FIXED_CONDITIONS_FINISHED !<Is .TRUE. if the fixed conditions for the problem has finished being created, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_BOUNDARY_CONDITIONS(:) !<GLOBAL_BOUNDARY_CONDITIONS(dof_idx). The global boundary condition for the dof_idx'th dof of the dependent field. \see PROBLEM_ROUTINES_FixedConditions
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the distributed vector containing the boundary conditions for the domain for this process.
  END TYPE PROBLEM_FIXED_CONDITIONS_TYPE

  !>Contains information on the dependent variables for the problem.
  TYPE PROBLEM_DEPENDENT_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem.
    LOGICAL :: DEPENDENT_FINISHED !<Is .TRUE. if the dependent variables for the problem has finished being created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<A pointer to the dependent field for the problem.
  END TYPE PROBLEM_DEPENDENT_TYPE

  !>Contains information on the source for the problem.
  TYPE PROBLEM_SOURCE_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem. 
    LOGICAL :: SOURCE_FINISHED !<Is .TRUE. if the source for the problem has finished being created, .FALSE. if not.
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD !<A pointer to the source field for the problem if one is defined. If no source is defined the pointer is NULL.
  END TYPE PROBLEM_SOURCE_TYPE

  !>Contains information on the analytic setup for the problem.
  TYPE PROBLEM_ANALYTIC_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem.
    LOGICAL :: ANALYTIC_FINISHED !<Is .TRUE. if the analytic setup for the problem has finished being created, .FALSE. if not.
  END TYPE PROBLEM_ANALYTIC_TYPE

  !>Contains information on the interpolation for the problem solution
  TYPE PROBLEM_INTERPOLATION_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field for the problem.
    TYPE(FIELD_TYPE), POINTER :: FIBRE_FIELD !<A pointer to the fibre field for the problem (if one is defined).
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<A pointer to the dependent field for the problem 
    TYPE(FIELD_TYPE), POINTER :: MATERIAL_FIELD !<A pointer to the material field for the problem (if one is defined).
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD !<A pointer to the source field for the problem (if one is defined).
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERP_PARAMETERS !<A pointer to the geometric interpolation parameters for the problem.
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FIBRE_INTERP_PARAMETERS !<A pointer to the fibre interpolation parameters for the problem (if a fibre field is defined). 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DEPENDENT_INTERP_PARAMETERS !<A pointer to the dependent interpolation parameters for the problem. 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: MATERIAL_INTERP_PARAMETERS !<A pointer to the material interpolation parameters for the problem (if a material field is defined). 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: SOURCE_INTERP_PARAMETERS !<A pointer to the source interpolation parameters for the problem (if a source field is defined). 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERP_POINT !<A pointer to the geometric interpolated point information for the problem. 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FIBRE_INTERP_POINT !<A pointer to the fibre interpolated point information for the problem (if a fibre field is defined). 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERP_POINT !<A pointer to the dependent interpolated point information for the problem. 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: MATERIAL_INTERP_POINT !<A pointer to the material interpolated point information for the problem (if a material field is defined). 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: SOURCE_INTERP_POINT !<A pointer to the source interpolated point information for the problem (if a source field is defined).
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERP_POINT_METRICS !<A pointer to the geometric interpolated point metrics information 
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: FIBRE_INTERP_POINT_METRICS !<A pointer to the fibre interpolated point metrics information 
  END TYPE PROBLEM_INTERPOLATION_TYPE

  !>Contains information on any data required for a linear solution
  TYPE PROBLEM_LINEAR_DATA_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution.
  END TYPE PROBLEM_LINEAR_DATA_TYPE

  !>Contains information on any data required for a non-linear solution
  TYPE PROBLEM_NONLINEAR_DATA_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution.
    INTEGER(INTG) :: NUMBER_OF_ITERATIONS
  END TYPE PROBLEM_NONLINEAR_DATA_TYPE

  !>Contains information on any data required for a time-dependent solution
  TYPE PROBLEM_TIME_DATA_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution.
  END TYPE PROBLEM_TIME_DATA_TYPE

  !>Contains information for an element matrix
  TYPE ELEMENT_MATRIX_TYPE
    INTEGER(INTG) :: GLOBAL_MATRIX_NUMBER
    INTEGER(INTG) :: NUMBER_OF_ROWS
    INTEGER(INTG) :: NUMBER_OF_COLUMNS
    INTEGER(INTG) :: MAX_NUMBER_OF_ROWS
    INTEGER(INTG) :: MAX_NUMBER_OF_COLUMNS
    INTEGER(INTG), ALLOCATABLE :: ROW_DOFS(:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_DOFS(:)
    REAL(DP), ALLOCATABLE :: MATRIX(:,:)
  END TYPE ELEMENT_MATRIX_TYPE

  !>Contains information for an element vector
  TYPE ELEMENT_VECTOR_TYPE
    INTEGER(INTG) :: NUMBER_OF_ROWS
    INTEGER(INTG) :: MAX_NUMBER_OF_ROWS
    INTEGER(INTG), ALLOCATABLE :: ROW_DOFS(:)
    REAL(DP), ALLOCATABLE :: VECTOR(:)
  END TYPE ELEMENT_VECTOR_TYPE
  
  !>Contains information on the mapping from a global matrix row/column to a solver matrix row/column.
  TYPE GLOBAL_TO_SOLVER_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_SOLUTION_DOFS !<Number of solver matrix row/column dofs the global matrix row/column dof is mapped to
    INTEGER(INTG), ALLOCATABLE :: SOLUTION_DOFS(:) !<SOLUTION_DOFS(i). Contains the i'th solution row/column dof that this global matrix row/column dof is mapped to.
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(i). Contains the i'th coupling coefficient
  END TYPE GLOBAL_TO_SOLVER_MAP_TYPE

  !>Contains information about a global matrix
  TYPE PROBLEM_GLOBAL_MATRIX_TYPE
    INTEGER(INTG) :: MATRIX_NUMBER !<The number of the global matrix
    INTEGER(INTG) :: VARIABLE_TYPE !<The variable type for this matrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE !<A pointer to the dependent field variable for this matrix
    INTEGER(INTG) :: STORAGE_TYPE !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: STRUCTURE_TYPE !<The structure (sparsity) type for this matrix
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if this global matrix is to be updated
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of rows in this global matrix
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in this global matrix
    INTEGER(INTG), ALLOCATABLE :: DOF_ROW_MAP(:) !<DOF_ROW_MAP(row_idx). The DOF that the row_idx'th row of this global matrix is mapped to.
    INTEGER(INTG), ALLOCATABLE :: DOF_COLUMN_MAP(:) !<DOF_COLUMN_MAP(column_idx). The DOF that the column_idx'th column of this global matrix is mapped to. 
    TYPE(GLOBAL_TO_SOLVER_MAP_TYPE), ALLOCATABLE :: SOLVER_ROW_MAP(:) !<SOLVER_ROW_MAP(row_idx). The mapping from the row_idx'th row of this global matrix to the solver matrix rows. 
    TYPE(GLOBAL_TO_SOLVER_MAP_TYPE), ALLOCATABLE :: SOLVER_COLUMN_MAP(:) !<SOLVER_COLUMN_MAP(column_idx). The mapping from the column_idx'th column of this global matrix to the solver matrix columns.
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the distributed global matrix data
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix for this glboal matrix
  END TYPE PROBLEM_GLOBAL_MATRIX_TYPE

  !>A buffer type to allow for an array of pointers to a GLOBAL_MATRIX_TYPE.
  TYPE PROBLEM_GLOBAL_MATRIX_PTR_TYPE    
    TYPE(PROBLEM_GLOBAL_MATRIX_TYPE), POINTER :: PTR !<The pointer to the global matrix
  END TYPE PROBLEM_GLOBAL_MATRIX_PTR_TYPE

  !>A type to temporarily hold (cache) the user modifiable values that are used to create a field.
  TYPE PROBLEM_GLOBAL_MATRICES_CREATE_VALUES_CACHE_TYPE
    INTEGER(INTG), ALLOCATABLE :: MATRIX_VARIABLE_TYPES(:) !<MATRIX_VARIABLE_TYPES(matrix_idx). The dependent variable type mapped to the matrix_idx'th global matrix
    INTEGER(INTG), ALLOCATABLE :: MATRIX_STORAGE_TYPE(:) !<MATRIX_STORAGE_TYPE(matrix_idx). The storage type (sparsity) for the matrix_idx'th global matrix.
    INTEGER(INTG), ALLOCATABLE :: MATRIX_STRUCTURE_TYPE(:) !<MATRIX_STRUCTURE_TYPE(matrix_idx). The structure type (sparsity) for the matrix_idx'th global matrix.
  END TYPE PROBLEM_GLOBAL_MATRICES_CREATE_VALUES_CACHE_TYPE
  
  !>Contains the information about the mapping of a DOF to global matrices
  TYPE DOF_TO_GLOBAL_MAPPING_TYPE
    INTEGER(INTG) :: ROW_MAPPING !<The global row number the DOF is mapped to.
    INTEGER(INTG) :: NUMBER_OF_MATRICES !<The number of global matrices the DOF is mapped to. Note that if the number of matrices is zero then the DOF is mapped to the global rhs vector. 
    INTEGER(INTG), ALLOCATABLE :: MATRIX_MAPPINGS(:) !<MATRIX_MAPPINGS(i). The i'th global matrix the DOF is mapped to
    INTEGER(INTG), ALLOCATABLE :: COLUMN_MAPPINGS(:) !<COLUMN_MAPPINGS(i). The global column number in the i'th global matrix that the DOF is mapped to.  
  END TYPE DOF_TO_GLOBAL_MAPPING_TYPE

  !>Contains information on the global matrices and rhs vector
  TYPE PROBLEM_GLOBAL_MATRICES_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer back to the problem soluion
    LOGICAL :: GLOBAL_MATRICES_FINISHED !<Is .TRUE. if the global matrices have finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of rows in the distributed global matrices
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ROWS !<The total number of rows in the distributed global matrices
    INTEGER(INTG) :: NUMBER_OF_MATRICES !<The number of global matrices defined for the problem.
    INTEGER(INTG) :: RHS_VARIABLE_TYPE !<The dependent variable type mapped to the global RHS vector. If the problem has no global RHS vector then the RHS_VARIABLE_TYPE is 0.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: RHS_VARIABLE !<A pointer to the dependent variable  mapped to the global RHS vector. If the problem has no global RHS vector then the RHS_VARIABLE is NULL.
    TYPE(PROBLEM_GLOBAL_MATRIX_PTR_TYPE), ALLOCATABLE :: MATRIX_TYPE_MAP(:) !<The matrix type map. 
    TYPE(PROBLEM_GLOBAL_MATRIX_TYPE), ALLOCATABLE :: MATRICES(:) !<MATRICES(matrix_idx) contains the information on the matrix_idx'th global matrix
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the global rhs vector is to be updated
    TYPE(DOF_TO_GLOBAL_MAPPING_TYPE), ALLOCATABLE :: DOF_TO_GLOBAL_MAPPING(:) !<DOF_TO_GLOBAL_MAPPING(ny). The mappings to the global matrices for the ny'th dependent DOF.
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ROW_TO_DOF_MAP(:) !<GLOBAL_ROW_TO_DOF_MAP(ny). The global DOF corresponding to the ny'th row of the global rhs vector. 
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the distributed global rhs vector data
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element rhs information
    TYPE(PROBLEM_GLOBAL_MATRICES_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<The create values cache for the global matrices
  END TYPE PROBLEM_GLOBAL_MATRICES_TYPE

  !>Contains the information on what global matrix rows contributed to this solver matrix row.
  TYPE SOLVER_TO_GLOBAL_MAP_TYPE
    INTEGER(INTG) :: NUMBER_OF_MATRICES !<The number of global matrices that contributed to this solution row.
    
    INTEGER(INTG), ALLOCATABLE :: COLUMNS(:) !<COLUMNS(matrix_idx) is the column number in the matrix_idx'th global matrix that this solution column/variable maps to
    REAL(DP), ALLOCATABLE :: COUPLING_COEFFICIENTS(:) !<COUPLING_COEFFICIENTS(matrix_idx) is the coupling coefficient for the matrix_idth'th global 
  END TYPE SOLVER_TO_GLOBAL_MAP_TYPE

  !>Contains information on the solver matrices and rhs vector
  TYPE PROBLEM_SOLVER_MATRICES_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution
    INTEGER(INTG) :: NUMBER_OF_ROWS !<The number of rows in the distributed solution matrix for this computational node
    INTEGER(INTG) :: TOTAL_NUMBER_OF_ROWS !<The total number of rows in the distributed solution matrix
    INTEGER(INTG) :: NUMBER_OF_COLUMNS !<The number of columns in the distributed solution matrix
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if the solution matrix is to be updated
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the solution vector is to be updated
    TYPE(SOLVER_TO_GLOBAL_MAP_TYPE), ALLOCATABLE :: SOLVER_TO_GLOBAL_MAP(:) !<SOLVER_TO_GLOBAL_MAP(no) is the solver to global mappings for the no'th solver column.
  END TYPE PROBLEM_SOLVER_MATRICES_TYPE

  !>Contains information on the type of solver to be used
  TYPE PROBLEM_SOLVER_TYPE
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution
    LOGICAL :: SOLVER_FINISHED !<Is .TRUE. if the problem solver has finished being created, .FALSE. if not.
  END TYPE PROBLEM_SOLVER_TYPE

  !>Contains information regarding the solution of a problem
  TYPE PROBLEM_SOLUTION_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    LOGICAL :: SOLUTION_FINISHED !<Is .TRUE. if the problem solution has finished being created, .FALSE. if not.
    INTEGER(INTG) :: OUTPUT_TYPE !<The output type for the problem solution \see PROBLEM_ROUTINES_SolutionOutputTypes,PROBLEM_ROUTINES
    TYPE(PROBLEM_INTERPOLATION_TYPE), POINTER :: INTERPOLATION !<A pointer to the interpolation information used in the problem solution
    TYPE(PROBLEM_LINEAR_DATA_TYPE), POINTER :: LINEAR_DATA !<A pointer to the data for linear problems.
    TYPE(PROBLEM_NONLINEAR_DATA_TYPE), POINTER :: NONLINEAR_DATA !<A pointer to the data for non-linear problems.
    TYPE(PROBLEM_TIME_DATA_TYPE), POINTER :: TIME_DATA !<A pointer to the data for non-static problems
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver used for the problem solution
    TYPE(PROBLEM_GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices and vectors used for the problem solution.
    TYPE(PROBLEM_SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices and vectors used for the problem solution.    
  END TYPE PROBLEM_SOLUTION_TYPE

  !>Contains information for a problem defined on a region.
  TYPE PROBLEM_TYPE
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the problem. The user number must be unique.
    INTEGER(INTG) :: GLOBAL_NUMBER !<The global number of the problem in the list of problems for a region.
    LOGICAL :: PROBLEM_FINISHED !<Is .TRUE. if the problem has finished being created, .FALSE. if not.
    TYPE(PROBLEMS_TYPE), POINTER :: PROBLEMS !<A pointer to the problems for this problem.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region.
    
    INTEGER(INTG) :: CLASS !<The problem specification class identifier
    INTEGER(INTG) :: TYPE !<The problem specification type identifier
    INTEGER(INTG) :: SUBTYPE !<The problem specification subtype identifier
    
    INTEGER(INTG) :: LINEARITY !<The problem linearity type \see PROBLEM_ROUTINES_LinearityTypes
    INTEGER(INTG) :: TIME_TYPE !<The problem time dependence type \see PROBLEM_ROUTINES_TimeDepedenceTypes
    
    INTEGER(INTG) :: SOLUTION_METHOD !<The solution method for the problem \see PROBLEM_ROUTINES_SolutionMethods 
    
    TYPE(PROBLEM_GEOMETRY_TYPE) :: GEOMETRY !<The geometry information for the problem.
    TYPE(PROBLEM_MATERIALS_TYPE), POINTER :: MATERIALS !<A pointer to the materials information for the problem.
    TYPE(PROBLEM_SOURCE_TYPE), POINTER :: SOURCE !<A pointer to the source information for the problem.
    TYPE(PROBLEM_DEPENDENT_TYPE) :: DEPENDENT !<The depedent variable information for the problem.
    TYPE(PROBLEM_ANALYTIC_TYPE), POINTER :: ANALYTIC !<A pointer to the analytic setup information for the problem.
    TYPE(PROBLEM_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS !<A pointer to the fixed condition information for the problem. \todo Change name to BOUNDARY_CONDITIONS???
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution information for the problem.
  END TYPE PROBLEM_TYPE
  
  !>A buffer type to allow for an array of pointers to a PROBLEM_TYPE.
  TYPE PROBLEM_PTR_TYPE
    TYPE(PROBLEM_TYPE), POINTER :: PTR !<The pointer to the problem.
  END TYPE PROBLEM_PTR_TYPE
       
  !>Contains information on the problems defined on a region.
  TYPE PROBLEMS_TYPE
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region.
    INTEGER(INTG) :: NUMBER_OF_PROBLEMS !<The number of problems defined on the region.
    TYPE(PROBLEM_PTR_TYPE), POINTER :: PROBLEMS(:) !<The array of pointers to the problems.
  END TYPE PROBLEMS_TYPE
  
  !>A buffer type to allow for an array of pointers to a REGION_TYPE.
  TYPE REGION_PTR_TYPE
    TYPE(REGION_TYPE), POINTER :: PTR !<The pointer to the region.
  END TYPE REGION_PTR_TYPE
     
  !>Contains information for a region.
  TYPE REGION_TYPE 
    INTEGER(INTG) :: USER_NUMBER !<The user defined identifier for the region. The user number must be unique.
    LOGICAL :: REGION_FINISHED !<Is .TRUE. if the region has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: LABEL !<A user defined label for the region.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system used by the region.
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes defined on the region.
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes defined on the region.
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointer to the fields defined on the region.
    TYPE(PROBLEMS_TYPE), POINTER :: PROBLEMS !<A pointer to the problems defined on the region.
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region for the region. If the region has no parent region then it is the global (world) region and PARENT_REGION is NULL.
    INTEGER(INTG) :: NUMBER_OF_SUB_REGIONS !<The number of sub-regions defined for the region.
    TYPE(REGION_PTR_TYPE), POINTER :: SUB_REGIONS(:) !<An array of pointers to the sub-regions defined on the region.
  END TYPE REGION_TYPE

END MODULE TYPES
