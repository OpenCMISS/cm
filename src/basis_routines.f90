!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all basis function routines.
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

!> This module contains all basis function routines.
MODULE BASIS_ROUTINES

  USE BASE_ROUTINES
  USE CONSTANTS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BASIS_ROUTINES_BasisTypes BASIS_ROUTINES::BasisTypes
  !> \brief Basis definition type parameters
  !> \todo Combine simplex and serendipity elements???
  !> \see BASIS_ROUTINES,OPENCMISS_BasisTypes
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LAGRANGE_HERMITE_TP_TYPE=1 !<Lagrange-Hermite tensor product basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SIMPLEX_TYPE=2 !<Simplex basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SERENDIPITY_TYPE=3 !<Serendipity basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_AUXILLIARY_TYPE=4 !<Auxillary basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_B_SPLINE_TP_TYPE=5 !<B-spline basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE=6 !<Fourier-Lagrange tensor product basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_EXTENDED_LAGRANGE_TP_TYPE=7 !< Extendend Lagrange tensor product basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_InterpolationSpecifications BASIS_ROUTINES::InterpolationSpecifications
  !> \brief Interpolation specification parameters
  !> \see BASIS_ROUTINES,OPENCMISS_InterpolationSpecifications
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_LAGRANGE_INTERPOLATION=1 !<Linear Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_LAGRANGE_INTERPOLATION=2 !<Quadratic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_LAGRANGE_INTERPOLATION=3 !<Cubic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_HERMITE_INTERPOLATION=4 !<Cubic Hermite interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC1_HERMITE_INTERPOLATION=5 !<Quadratic Hermite (no derivative at xi=0) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC2_HERMITE_INTERPOLATION=6 !<Quadratic Hermite (no derivative at xi=1) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_SIMPLEX_INTERPOLATION=7 !<Linear Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_SIMPLEX_INTERPOLATION=8 !<Quadratic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_SIMPLEX_INTERPOLATION=9 !<Cubic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_InterpolationTypes BASIS_ROUTINES::InterpolationTypes
  !> \brief Interpolation type parameters for a Xi direction
  !> \see BASIS_ROUTINES
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LAGRANGE_INTERPOLATION=1 !<Lagrange interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_HERMITE_INTERPOLATION=2 !<Hermite interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SIMPLEX_INTERPOLATION=3 !<Simplex interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SERENDIPITY_INTERPOLATION=4 !<Serendipity interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_TRANSITION_INTERPOLATION=5 !<Transition interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SINGULAR_INTERPOLATION=6 !<Singular interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_FOURIER_INTERPOLATION=7 !<Fourier interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  !>@}
  
  !> \addtogroup BASIS_ROUTINES_InterpolationOrder BASIS_ROUTINES::InterpolationOrder
  !> \brief Interpolation order for a Xi direction
  !> \see BASIS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_INTERPOLATION_ORDER=1 !<Linear interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_INTERPOLATION_ORDER=2 !<Quadratic interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_INTERPOLATION_ORDER=3 !<Cubic interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC1_INTERPOLATION_ORDER=4 !<Quadratic (no derivative at xi=0) interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC2_INTERPOLATION_ORDER=5 !<Quadratic (no derivative at xi=1) interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_QuadratureSchemes BASIS_ROUTINES::QuadratureSchemes
  !> \brief Quadrature scheme parameters
  !> \see BASIS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES=4 !<The number of currently defined quadrature schemes \see BASIS_ROUTINES_QuadratureSchemes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_DEFAULT_QUADRATURE_SCHEME=1 !<Identifier for the default quadrature scheme \see BASIS_ROUTINES_QuadratureSchemes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_LOW_QUADRATURE_SCHEME=2 !<Identifier for a low order quadrature scheme \see BASIS_ROUTINES_QuadratureSchemes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_MID_QUADRATURE_SCHEME=3 !<Identifier for a mid order quadrature scheme \see BASIS_ROUTINES_QuadratureSchemes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_HIGH_QUADRATURE_SCHEME=4 !<Identifier for a high order quadrature scheme \see BASIS_ROUTINES_QuadratureSchemes,BASIS_ROUTINES
  !>@}
  
  !> \addtogroup BASIS_ROUTINES_QuadratureTypes BASIS_ROUTINES::QuadratureTypes
  !> \brief Quadrature type parameters
  !> \see BASIS_ROUTINES,OPENCMISS_QuadratureTypes
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_LEGENDRE_QUADRATURE=1 !<Gauss-Legendre quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_LAGUERRE_QUADRATURE=2 !<Gauss-Laguerre quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GUASS_HERMITE_QUADRATURE=3 !<Gauss-Hermite quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE=4 !<Adaptive Gauss-Legendre quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_SIMPLEX_QUADRATURE=5 !<Gauss-Legendre for Simplex elements quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_XiCollapse BASIS_ROUTINES::XiCollapse
  !> \brief Xi collapse parameters
  !> \see BASIS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_XI_COLLAPSED=1 !<The Xi direction is collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_COLLAPSED_AT_XI0=2 !<The Xi direction at the xi=0 end of this Xi direction is collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_COLLAPSED_AT_XI1=3 !<The Xi direction at the xi=1 end of this Xi direction is collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_NOT_COLLAPSED=4 !<The Xi direction is not collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  !>@}
  
  !!Module types
  ! 
  !!>Contains information on the defined basis functions
  !TYPE BASIS_FUNCTIONS_TYPE
  !  PRIVATE
  !  INTEGER(INTG) :: NUMBER_BASIS_FUNCTIONS !<The number of basis functions definegd
  !  TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:) !<The array of pointers to the defined basis functions
  !END TYPE BASIS_FUNCTIONS_TYPE
  
  !Module variables

  TYPE(BASIS_FUNCTIONS_TYPE) :: BASIS_FUNCTIONS !<The tree of defined basis functions
  
  !Interfaces

  !>Evaluates the appropriate partial derivative index for the specificied basis function at a Xi location \see BASIS_ROUTINES
  INTERFACE BASIS_EVALUATE_XI
    MODULE PROCEDURE BASIS_EVALUATE_XI_DP
  END INTERFACE !BASIS_EVALUATE_XI

  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at a Gauss point \see BASIS_ROUTINES
  INTERFACE BASIS_INTERPOLATE_GAUSS
    MODULE PROCEDURE BASIS_INTERPOLATE_GAUSS_DP
  END INTERFACE !BASIS_INTERPOLATE_GAUSS
  
  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at Xi location \see BASIS_ROUTINES
  INTERFACE BASIS_INTERPOLATE_XI
    MODULE PROCEDURE BASIS_INTERPOLATE_XI_DP
  END INTERFACE !BASIS_INTERPOLATE_XI

  !>Interpolates the requested partial derivative index(ices) of the element parameters for basis function at a face Gauss point \see BASIS_ROUTINES
  INTERFACE BASIS_INTERPOLATE_LOCAL_FACE_GAUSS
    MODULE PROCEDURE BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP
  END INTERFACE !BASIS_INTERPOLATE_LOCAL_FACE_GAUSS

  !>Sets/changes the interpolation type in each Xi direction for a basis
  INTERFACE BASIS_INTERPOLATION_XI_SET
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_SET_PTR
  END INTERFACE !BASIS_INTERPOLATION_XI_SET

  !>Sets/changes the number of Xi directions for a basis
  INTERFACE BASIS_NUMBER_OF_XI_SET
    MODULE PROCEDURE BASIS_NUMBER_OF_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_NUMBER_OF_XI_SET_PTR
  END INTERFACE !BASIS_NUMBER_OF_XI_SET

  !>Sets/changes the type for a basis.
  INTERFACE BASIS_TYPE_SET
    MODULE PROCEDURE BASIS_TYPE_SET_NUMBER
    MODULE PROCEDURE BASIS_TYPE_SET_PTR
  END INTERFACE !BASIS_TYPE_SET

  !>Sets/changes the collapsed Xi flags for a basis.
  INTERFACE BASIS_COLLAPSED_XI_SET
    MODULE PROCEDURE BASIS_COLLAPSED_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_COLLAPSED_XI_SET_PTR
  END INTERFACE !BASIS_COLLAPSED_XI_SET

  !>Sets/changes the number of gauss in each Xi direction for a basis quadrature.
  INTERFACE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET
    MODULE PROCEDURE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR
  END INTERFACE !BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET

  !>Sets/changes the order of a quadrature for a basis quadrature.
  INTERFACE BASIS_QUADRATURE_ORDER_SET
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_SET_PTR
  END INTERFACE !BASIS_QUADRATURE_ORDER_SET

  !>Sets/changes the quadrature type for a basis
  INTERFACE BASIS_QUADRATURE_TYPE_SET
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_SET_PTR
  END INTERFACE !BASIS_QUADRATURE_TYPE_SET
  
  INTERFACE SIMPLEX_LINEAR_EVALUATE
    MODULE PROCEDURE SIMPLEX_LINEAR_EVALUATE_DP
  END INTERFACE !SIMPLEX_LINEAR_EVALUATE

  INTERFACE SIMPLEX_QUADRATIC_EVALUATE
    MODULE PROCEDURE SIMPLEX_QUADRATIC_EVALUATE_DP
  END INTERFACE !SIMPLEX_QUADRATIC_EVALUATE

  INTERFACE SIMPLEX_CUBIC_EVALUATE
    MODULE PROCEDURE SIMPLEX_CUBIC_EVALUATE_DP
  END INTERFACE !SIMPLEX_CUBIC_EVALUATE

  !>Evaluates the Lagrange/Hermite/Fourier tensor product basis function for the given basis
  INTERFACE BASIS_LHTP_BASIS_EVALUATE
    MODULE PROCEDURE BASIS_LHTP_BASIS_EVALUATE_DP
  END INTERFACE !BASIS_LHTP_BASIS_EVALUATE

  PUBLIC BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE,BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
    & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE

  PUBLIC BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & BASIS_CUBIC_HERMITE_INTERPOLATION,BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION, &
    & BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,BASIS_CUBIC_SIMPLEX_INTERPOLATION

  PUBLIC BASIS_LINEAR_INTERPOLATION_ORDER,BASIS_QUADRATIC_INTERPOLATION_ORDER,BASIS_CUBIC_INTERPOLATION_ORDER, &
    & BASIS_QUADRATIC1_INTERPOLATION_ORDER,BASIS_QUADRATIC2_INTERPOLATION_ORDER
  
  PUBLIC BASIS_DEFAULT_QUADRATURE_SCHEME,BASIS_LOW_QUADRATURE_SCHEME,BASIS_MID_QUADRATURE_SCHEME,BASIS_HIGH_QUADRATURE_SCHEME
  
  PUBLIC BASIS_GAUSS_LEGENDRE_QUADRATURE,BASIS_GAUSS_LAGUERRE_QUADRATURE,BASIS_GUASS_HERMITE_QUADRATURE,&
    & BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE,BASIS_GAUSS_SIMPLEX_QUADRATURE

  PUBLIC BASIS_XI_COLLAPSED,BASIS_COLLAPSED_AT_XI0,BASIS_COLLAPSED_AT_XI1,BASIS_NOT_COLLAPSED, BASIS_FUNCTIONS

  PUBLIC BASIS_COLLAPSED_XI_SET
  
  PUBLIC BASIS_EVALUATE_XI
  
  PUBLIC BASIS_INTERPOLATE_GAUSS,BASIS_INTERPOLATE_XI,BASIS_INTERPOLATE_LOCAL_FACE_GAUSS

  PUBLIC BASIS_LOCAL_NODE_XI_CALCULATE

  PUBLIC BASIS_NUMBER_OF_LOCAL_NODES_GET

  PUBLIC BASIS_INTERPOLATION_XI_SET

  PUBLIC BASIS_NUMBER_OF_XI_SET

  PUBLIC BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET

  PUBLIC BASIS_QUADRATURE_DESTROY,BASIS_QUADRATURE_ORDER_SET,BASIS_QUADRATURE_TYPE_SET

  PUBLIC BASIS_TYPE_SET,BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET

  PUBLIC BASIS_CREATE_START,BASIS_CREATE_FINISH

  PUBLIC BASIS_DESTROY

  PUBLIC BASES_FINALISE,BASES_INITIALISE

  PUBLIC BASIS_USER_NUMBER_FIND
    
  PUBLIC BASIS_COLLAPSED_XI_GET,BASIS_INTERPOLATION_XI_GET,BASIS_NUMBER_OF_XI_GET,BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET, &
    & BASIS_QUADRATURE_ORDER_GET,BASIS_QUADRATURE_TYPE_GET,BASIS_TYPE_GET



CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalises the bases and deallocates all memory
  SUBROUTINE BASES_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASES_FINALISE",ERR,ERROR,*999)

    !Destroy any created basis functions
    DO WHILE(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS>0)
      CALL BASIS_DESTROY(BASIS_FUNCTIONS%BASES(1)%PTR,ERR,ERROR,*999)
    ENDDO !nb
    !Destroy basis functions and deallocated any memory allocated
    BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS=0
    IF(ASSOCIATED(BASIS_FUNCTIONS%BASES)) DEALLOCATE(BASIS_FUNCTIONS%BASES)
    
    CALL EXITS("BASES_FINALISE")
    RETURN
999 CALL ERRORS("BASES_FINALISE",ERR,ERROR)
    CALL EXITS("BASES_FINALISE")
    RETURN 1
    
  END SUBROUTINE BASES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the bases.
  SUBROUTINE BASES_INITIALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASES_INITIALISE",ERR,ERROR,*999)

    BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS=0    
    NULLIFY(BASIS_FUNCTIONS%BASES)
   
    CALL EXITS("BASES_INITIALISE")
    RETURN
999 CALL ERRORS("BASES_INITIALISE",ERR,ERROR)
    CALL EXITS("BASES_INITIALISE")
    RETURN 1
    
  END SUBROUTINE BASES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a new basis \see BASIS_ROUTINES::BASIS_CREATE_START,OPENCMISS::BasisCreateFinish
  SUBROUTINE BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni,nic,nn,nn1,nn2,nn3,nn4,ns,local_line_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      SELECT CASE(BASIS%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        CALL BASIS_LHTP_FAMILY_CREATE(BASIS,ERR,ERROR,*999)
      CASE(BASIS_SIMPLEX_TYPE)
        CALL BASIS_SIMPLEX_FAMILY_CREATE(BASIS,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Basis type "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))//" is invalid or not implemented"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      BASIS%BASIS_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Basis user number = ",BASIS%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Basis family number = ",BASIS%FAMILY_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Basis global number = ",BASIS%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Basis type = ",BASIS%TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Basis degenerate = ",BASIS%DEGENERATE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi directions = ",BASIS%NUMBER_OF_XI,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES,ERR,ERROR,*999)
      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI_COORDINATES,4,4,BASIS%INTERPOLATION_TYPE, &
        & '("  Interpolation type(nic):",4(X,I2))','(25X,4(X,I2))',ERR,ERROR,*999)      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI_COORDINATES,4,4,BASIS%INTERPOLATION_ORDER, &
        & '("  Interpolation order(nic):",4(X,I2))','(26X,4(X,I2))',ERR,ERROR,*999)      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI,3,3,BASIS%COLLAPSED_XI, &
        & '("  Collapsed xi(ni):",3(X,I2))','(26X,3(X,I2))',ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of partial derivatives = ",BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of nodes = ",BASIS%NUMBER_OF_NODES,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI_COORDINATES,4,4,BASIS%NUMBER_OF_NODES_XIC, &
        & '("  Number of nodes(nic):",4(X,I2))','(22X,4(X,I2))',ERR,ERROR,*999)      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_NODES,8,8,BASIS%NUMBER_OF_DERIVATIVES, &
        & '("  Number of derivatives(nn):",8(X,I2))','(28X,8(X,I2))',ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of element parameters = ",BASIS%NUMBER_OF_ELEMENT_PARAMETERS, &
        & ERR,ERROR,*999)
! CPB 23/07/07 Doxygen may or may not like this line for some reason????      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_NODES,8,8,BASIS%NODE_AT_COLLAPSE, &
        & '("  Node at collapse(nn):",8(X,L1))','(23X,8(X,L1))',ERR,ERROR,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Node position index:",ERR,ERROR,*999)
      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Xic = ",nic,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_NODES,16,16,BASIS%NODE_POSITION_INDEX(:,nic), &
          & '("      INDEX(nn)  :",16(X,I2))','(18X,16(X,I2))',ERR,ERROR,*999)        
      ENDDO !ni
      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse node position index:",ERR,ERROR,*999)
      SELECT CASE(BASIS%NUMBER_OF_XI_COORDINATES)
      CASE(1)
        DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Xic = 1, Node position index = ",nn1,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      INDEX = ",BASIS%NODE_POSITION_INDEX_INV(nn1,1,1,1), &
            & ERR,ERROR,*999)          
        ENDDO !nn1
      CASE(2)
        DO nn2=1,BASIS%NUMBER_OF_NODES_XIC(2)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Xic = 2, Node position index = ",nn2,ERR,ERROR,*999)
          DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xic = 1, Node position index = ",nn1,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        INDEX = ",BASIS%NODE_POSITION_INDEX_INV(nn1,nn2,1,1), &
              & ERR,ERROR,*999)
          ENDDO !nn1
        ENDDO !nn2
      CASE(3)
        DO nn3=1,BASIS%NUMBER_OF_NODES_XIC(3)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Xic = 3, Node position index = ",nn3,ERR,ERROR,*999)
          DO nn2=1,BASIS%NUMBER_OF_NODES_XIC(2)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xic = 2, Node position index = ",nn2,ERR,ERROR,*999)
            DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Xic = 1, Node position index = ",nn1,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          INDEX = ",BASIS%NODE_POSITION_INDEX_INV(nn1,nn2,nn3,1), &
                & ERR,ERROR,*999)
            ENDDO !nn1
          ENDDO !nn2
        ENDDO !nn3
      CASE(4)
        DO nn4=1,BASIS%NUMBER_OF_NODES_XIC(4)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Xic = 4, Node position index = ",nn4,ERR,ERROR,*999)
          DO nn3=1,BASIS%NUMBER_OF_NODES_XIC(3)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xic = 3, Node position index = ",nn3,ERR,ERROR,*999)
            DO nn2=1,BASIS%NUMBER_OF_NODES_XIC(2)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Xic = 2, Node position index = ",nn2,ERR,ERROR,*999)
              DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1) 
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Xic = 1, Node position index = ",nn1,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            INDEX = " &
                  & ,BASIS%NODE_POSITION_INDEX_INV(nn1,nn2,nn3,nn4),ERR,ERROR,*999)
              ENDDO !nn1
            ENDDO !nn2
          ENDDO !nn3
        ENDDO !nn4
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of xi coordinates",ERR,ERROR,*999)
      END SELECT
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative order index:",ERR,ERROR,*999)
      DO ni=1,BASIS%NUMBER_OF_XI
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Xi = ",ni,ERR,ERROR,*999)
        DO nn=1,BASIS%NUMBER_OF_NODES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node = ",nn,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & BASIS%DERIVATIVE_ORDER_INDEX(:,nn,ni),'("        INDEX(nk):",8(X,I2))','(18X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
      ENDDO !ni
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse derivative order index:",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Element parameter index:",ERR,ERROR,*999)
      DO nn=1,BASIS%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node = ",nn,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
          & BASIS%ELEMENT_PARAMETER_INDEX(:,nn),'("      INDEX(nk)  :",8(X,I2))','(18X,8(X,I2))',ERR,ERROR,*999)
      ENDDO !nn
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse element parameter index:",ERR,ERROR,*999)
      DO ns=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Element parameter = ",ns,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2, &
          & BASIS%ELEMENT_PARAMETER_INDEX_INV(:,ns),'("      INDEX(:)  :",2(X,I2))','(18X,2(X,I2))',ERR,ERROR,*999)
      ENDDO !ns
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Partial derivative index:",ERR,ERROR,*999)
      DO nn=1,BASIS%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node = ",nn,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
          & BASIS%PARTIAL_DERIVATIVE_INDEX(:,nn),'("      INDEX(nk)  :",8(X,I2))','(18X,8(X,I2))',ERR,ERROR,*999)
      ENDDO !nn
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local lines:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local lines = ",BASIS%NUMBER_OF_LOCAL_LINES,ERR,ERROR,*999)
      DO local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local line = ",local_line_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Local line xi direction = ", &
          & BASIS%LOCAL_LINE_XI_DIRECTION(local_line_idx),ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of nodes in local line = ", &
          & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(local_line_idx),ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(local_line_idx),4,4, &
          & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(:,local_line_idx),'("      Nodes in local line       :",4(X,I2))','(33X,4(X,I2))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(local_line_idx),4,4, &
          & BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(:,local_line_idx),'("      Derivatives in local line :",4(X,I2))', &
          & '(33X,4(X,I2))',ERR,ERROR,*999)
        IF(BASIS%NUMBER_OF_XI==2) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Local line xi normal = ", &
            & BASIS%LOCAL_XI_NORMAL(local_line_idx),ERR,ERROR,*999)
        ENDIF
      ENDDO !ni
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of sub-bases = ",BASIS%NUMBER_OF_SUB_BASES,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("BASIS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("BASIS_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE BASIS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a new basis 
  !>The default values of the BASIS attributes are:
  !>- TYPE: 1 (BASIS_LAGRANGE_HERMITE_TP_TYPE)
  !>- NUMBER_OF_XI: 3
  !>- INTERPOLATION_XI: (1,1,1) (BASIS_LINEAR_LAGRANGE_INTERPOLATIONs)
  !>- INTERPOLATION_TYPE: (1,1,1) (BASIS_LAGRANGE_INTERPOLATIONs)
  !>- INTERPOLATION_ORDER: (1,1,1) (BASIS_LINEAR_INTERPOLATION_ORDERs)
  !>- DEGENERATE: false
  !>- COLLAPSED_XI: (4,4,4) (BASIS_NOT_COLLAPSEDs)
  !>- QUADRATURE: 
  !>  - TYPE: 1 (BASIS_LAGRANGE_HERMITE_TP_TYPE)
  !>  - NUMBER_OF_GAUSS_XI: (2,2,2)
  !>  - GAUSS_ORDER: 0 
  SUBROUTINE BASIS_CREATE_START(USER_NUMBER,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to start the creation of
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the created basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nb
    TYPE(BASIS_TYPE), POINTER :: NEW_BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: NEW_BASES(:)

    NULLIFY(NEW_BASIS)
    NULLIFY(NEW_BASES)

    CALL ENTERS("BASIS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      CALL FLAG_ERROR("Basis is already associated",ERR,ERROR,*999)
    ELSE
      !See if basis number has already been created
      CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
      IF(ASSOCIATED(BASIS)) THEN
        CALL FLAG_ERROR("Basis number is already defined",ERR,ERROR,*999)
      ELSE
        !Allocate new basis function and add it to the basis functions
        ALLOCATE(NEW_BASES(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NEW_BASES",ERR,ERROR,*999)
        ALLOCATE(NEW_BASIS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NEW_BASIS",ERR,ERROR,*999)
        CALL BASIS_INITIALISE(NEW_BASIS,ERR,ERROR,*999)
        DO nb=1,BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS
          NEW_BASES(nb)%PTR=>BASIS_FUNCTIONS%BASES(nb)%PTR
        ENDDO !nb
        BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS=BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS+1
        NEW_BASES(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS)%PTR=>NEW_BASIS
        IF(ASSOCIATED(BASIS_FUNCTIONS%BASES)) DEALLOCATE(BASIS_FUNCTIONS%BASES)
        BASIS_FUNCTIONS%BASES=>NEW_BASES
        !Initialise the new basis function pointers
        NEW_BASIS%NUMBER_OF_SUB_BASES=0
        NULLIFY(NEW_BASIS%SUB_BASES)
        NULLIFY(NEW_BASIS%PARENT_BASIS)
        !Set the basis parameters
        NEW_BASIS%BASIS_FINISHED=.FALSE.
        NEW_BASIS%USER_NUMBER=USER_NUMBER
        NEW_BASIS%FAMILY_NUMBER=0
        NEW_BASIS%GLOBAL_NUMBER=BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS
        !Set the default basis parameters
        NEW_BASIS%TYPE=BASIS_LAGRANGE_HERMITE_TP_TYPE
        NEW_BASIS%NUMBER_OF_XI=3
        ALLOCATE(NEW_BASIS%INTERPOLATION_XI(3),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis interpolation xi",ERR,ERROR,*999)
        NEW_BASIS%INTERPOLATION_XI=(/BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
          & BASIS_LINEAR_LAGRANGE_INTERPOLATION/)
        ALLOCATE(NEW_BASIS%COLLAPSED_XI(3),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis collapsed xi",ERR,ERROR,*999)
        NEW_BASIS%COLLAPSED_XI=BASIS_NOT_COLLAPSED
        !Initialise the basis quadrature
        NULLIFY(NEW_BASIS%QUADRATURE%BASIS)
        CALL BASIS_QUADRATURE_INITIALISE(NEW_BASIS,ERR,ERROR,*999)
        
        BASIS=>NEW_BASIS
      ENDIF
    ENDIF
    
    CALL EXITS("BASIS_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_BASIS)) CALL BASIS_DESTROY(NEW_BASIS,ERR,ERROR,*998)
998 IF(ASSOCIATED(NEW_BASES)) DEALLOCATE(NEW_BASES)
    NULLIFY(BASIS)
    CALL ERRORS("BASIS_CREATE_START",ERR,ERROR)
    CALL EXITS("BASIS_CREATE_START")
    RETURN 1
    
  END SUBROUTINE BASIS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a basis identified by its basis user number \see BASIS_ROUTINES::BASIS_DESTROY_FAMILY,OPENCMISS::CMISSBasisDestroy
  RECURSIVE SUBROUTINE BASIS_DESTROY_NUMBER(USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASIS_DESTROY_NUMBER",ERR,ERROR,*999)

    CALL BASIS_FAMILY_DESTROY(USER_NUMBER,0,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_DESTROY_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_DESTROY_NUMBER")
    RETURN 1
    
  END SUBROUTINE BASIS_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a basis. \see BASIS_ROUTINES::BASIS_DESTROY_FAMILY,OPENCMISS::CMISSBasisDestroy
  RECURSIVE SUBROUTINE BASIS_DESTROY(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: USER_NUMBER
        
    CALL ENTERS("BASIS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      USER_NUMBER=BASIS%USER_NUMBER
      CALL BASIS_FAMILY_DESTROY(USER_NUMBER,0,ERR,ERROR,*999)
      !NULLIFY(BASIS)
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_DESTROY")
    RETURN
999 CALL ERRORS("BASIS_DESTROY",ERR,ERROR)
    CALL EXITS("BASIS_DESTROY")
    RETURN 1
    
  END SUBROUTINE BASIS_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the appropriate partial derivative index at position XI for the basis for double precision arguments.
  !>Note for simplex basis functions the XI coordinates should exclude the last area coordinate.
  FUNCTION BASIS_EVALUATE_XI_DP(BASIS,ELEMENT_PARAMETER_INDEX,PARTIAL_DERIV_INDEX,XI,ERR,ERROR)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: ELEMENT_PARAMETER_INDEX !<The element parameter index to evaluate i.e., the local basis index within the element basis.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to evaluate \see CONSTANTS_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: XI(:) !<The Xi position to evaluate the basis function at
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_EVALUATE_XI_DP
    !Local Variables
    INTEGER(INTG) :: nn,nk
    REAL(DP) :: XIL(SIZE(XI,1)+1)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BASIS_EVALUATE_XI_DP",ERR,ERROR,*999)
    
    BASIS_EVALUATE_XI_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(ELEMENT_PARAMETER_INDEX>0.AND.ELEMENT_PARAMETER_INDEX<=BASIS%NUMBER_OF_ELEMENT_PARAMETERS) THEN
        SELECT CASE(BASIS%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          nn=BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ELEMENT_PARAMETER_INDEX)
          nk=BASIS%ELEMENT_PARAMETER_INDEX_INV(2,ELEMENT_PARAMETER_INDEX)
          BASIS_EVALUATE_XI_DP=BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,PARTIAL_DERIV_INDEX,XI,ERR,ERROR)
          IF(ERR/=0) GOTO 999
        CASE(BASIS_SIMPLEX_TYPE)
          !Create the area coordinates from the xi coordinates
          XIL(1:SIZE(XI,1))=1.0_DP-XI
          XIL(SIZE(XI,1)+1)=SUM(XI)-(SIZE(XI,1)-1.0_DP)
          nn=BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ELEMENT_PARAMETER_INDEX)
          BASIS_EVALUATE_XI_DP=BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,PARTIAL_DERIV_INDEX,XIL,ERR,ERROR)
          IF(ERR/=0) GOTO 999
        CASE DEFAULT
          LOCAL_ERROR="Basis type "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))//" is invalid or not implemented."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The specified element parameter index of "// &
          & TRIM(NUMBER_TO_VSTRING(ELEMENT_PARAMETER_INDEX,"*",ERR,ERROR))// &
          & " is invalid. The index must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_EVALUATE_XI_DP")
    RETURN
999 CALL ERRORS("BASIS_EVALUATE_XI_DP",ERR,ERROR)
    CALL EXITS("BASIS_EVALUATE_XI_DP")
    RETURN
    
  END FUNCTION BASIS_EVALUATE_XI_DP
  
  !
  !================================================================================================================================
  !
  
  !>Destroys a basis identified by its basis user number and family number. Called from the library visible routine BASIS_DESTROY
  !> \see BASIS_ROUTINES::BASIS_DESTROY
  RECURSIVE SUBROUTINE BASIS_FAMILY_DESTROY(USER_NUMBER,FAMILY_NUMBER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to destroy
    INTEGER(INTG), INTENT(IN) :: FAMILY_NUMBER !<The family number of the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: count,nb
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: NEW_SUB_BASES(:)
    
    CALL ENTERS("BASIS_FAMILY_DESTROY",ERR,ERROR,*999)

    NULLIFY(BASIS)
    CALL BASIS_FAMILY_NUMBER_FIND(USER_NUMBER,FAMILY_NUMBER,BASIS,ERR,ERROR,*999)
    IF(ASSOCIATED(BASIS)) THEN

!!NOTE: We have to find a pointer to the basis to destroy within this routine rather than passing in a pointer to a
!!DESTROY_BASIS_PTR type routine because we need to change BASIS%SUB_BASES of the PARENT basis and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy BASIS pointer argument was associated with the SUB_BASES(x)%PTR actual
!!argument.
      
      IF(BASIS%NUMBER_OF_SUB_BASES==0) THEN
        !No more sub-bases so delete this instance
        IF(ASSOCIATED(BASIS%PARENT_BASIS)) THEN
          !Sub-basis function - delete this instance from the PARENT_BASIS
          NULLIFY(NEW_SUB_BASES)
          IF(BASIS%PARENT_BASIS%NUMBER_OF_SUB_BASES>1) THEN
            !If the parent basis has more than one sub basis then remove this instance from its sub-bases list
            ALLOCATE(NEW_SUB_BASES(BASIS%PARENT_BASIS%NUMBER_OF_SUB_BASES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-bases",ERR,ERROR,*999)
            count=0
            DO nb=1,BASIS%PARENT_BASIS%NUMBER_OF_SUB_BASES
              IF(BASIS%PARENT_BASIS%SUB_BASES(nb)%PTR%USER_NUMBER==BASIS%USER_NUMBER.AND. &
                & BASIS%PARENT_BASIS%SUB_BASES(nb)%PTR%FAMILY_NUMBER/=BASIS%FAMILY_NUMBER) THEN
                count=count+1
                NEW_SUB_BASES(count)%PTR=>BASIS%PARENT_BASIS%SUB_BASES(nb)%PTR
              ENDIF
            ENDDO
            IF(ASSOCIATED(BASIS%PARENT_BASIS%SUB_BASES)) DEALLOCATE(BASIS%PARENT_BASIS%SUB_BASES)
          ENDIF
          BASIS%PARENT_BASIS%NUMBER_OF_SUB_BASES=BASIS%PARENT_BASIS%NUMBER_OF_SUB_BASES-1
          BASIS%PARENT_BASIS%SUB_BASES=>NEW_SUB_BASES
        ELSE
          !Master basis function - delete this instance from BASIS_FUNCTIONS
          NULLIFY(NEW_SUB_BASES)
          IF(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS>1) THEN
            !If there is more than one basis defined then remove this instance from the basis functions
            ALLOCATE(NEW_SUB_BASES(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-bases",ERR,ERROR,*999)
            count=0
            DO nb=1,BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS
              IF(BASIS_FUNCTIONS%BASES(nb)%PTR%USER_NUMBER/=BASIS%USER_NUMBER.AND. &
                & BASIS_FUNCTIONS%BASES(nb)%PTR%FAMILY_NUMBER==0) THEN
                count=count+1
                NEW_SUB_BASES(count)%PTR=>BASIS_FUNCTIONS%BASES(nb)%PTR
              ENDIF
            ENDDO
            IF(ASSOCIATED(BASIS_FUNCTIONS%BASES)) DEALLOCATE(BASIS_FUNCTIONS%BASES)
          ENDIF
          BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS=BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS-1
          BASIS_FUNCTIONS%BASES=>NEW_SUB_BASES
        ENDIF

        CALL BASIS_FINALISE(BASIS,ERR,ERROR,*999)
         
      ELSE
        !Recursively delete sub-bases first
        DO WHILE(BASIS%NUMBER_OF_SUB_BASES>0)
          CALL BASIS_FAMILY_DESTROY(BASIS%SUB_BASES(1)%PTR%USER_NUMBER,BASIS%SUB_BASES(1)%PTR%FAMILY_NUMBER,ERR,ERROR,*999)
        ENDDO
        !Now delete this instance
        CALL BASIS_FAMILY_DESTROY(USER_NUMBER,FAMILY_NUMBER,ERR,ERROR,*999)
      ENDIF

    ELSE
      CALL FLAG_ERROR("Basis user number does not exist",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_FAMILY_DESTROY")
    RETURN
999 CALL ERRORS("BASIS_FAMILY_DESTROY",ERR,ERROR)
    CALL EXITS("BASIS_FAMILY_DESTROY")
    RETURN 1
  END SUBROUTINE BASIS_FAMILY_DESTROY

  !
  !================================================================================================================================
  !

  !>Finds and returns in BASIS a pointer to the basis with the given USER_NUMBER and FAMILY_NUMBER. If no basis with that
  !>number and family number exists then BASIS is returned nullified \see BASIS_ROUTINES::BASIS_USER_NUMBER_FIND
  RECURSIVE SUBROUTINE BASIS_FAMILY_NUMBER_FIND(USER_NUMBER,FAMILY_NUMBER,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to find
    INTEGER(INTG), INTENT(IN) :: FAMILY_NUMBER !<The family number of the basis to find
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On exit, A pointer to the basis. If no basis with the specified user and family numbers can be found the pointer is not associated.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nb,nsb
    TYPE(BASIS_TYPE), POINTER :: SUB_BASIS

    CALL ENTERS("BASIS_FAMILY_NUMBER_FIND",ERR,ERROR,*999)
    
    NULLIFY(BASIS)      
    nb=1
    DO WHILE(nb<=BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS.AND..NOT.ASSOCIATED(BASIS))
      IF(BASIS_FUNCTIONS%BASES(nb)%PTR%USER_NUMBER==USER_NUMBER) THEN
        IF(FAMILY_NUMBER==0) THEN
          BASIS=>BASIS_FUNCTIONS%BASES(nb)%PTR
        ELSE
!!TODO: \todo This only works for one level of sub-bases at the moment
          nsb=1
          DO WHILE(nsb<=BASIS_FUNCTIONS%BASES(nb)%PTR%NUMBER_OF_SUB_BASES.AND..NOT.ASSOCIATED(BASIS))
            SUB_BASIS=>BASIS_FUNCTIONS%BASES(nb)%PTR%SUB_BASES(nsb)%PTR
            IF(SUB_BASIS%FAMILY_NUMBER==FAMILY_NUMBER) THEN
              BASIS=>SUB_BASIS
            ELSE
              nsb=nsb+1
            ENDIF
          ENDDO
        ENDIF
      ELSE
        nb=nb+1
      ENDIF
    END DO
  
    CALL EXITS("BASIS_FAMILY_NUMBER_FIND")
    RETURN
999 CALL ERRORS("BASIS_FAMILY_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("BASIS_FAMILY_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE BASIS_FAMILY_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises a basis and deallocates all memory.
  SUBROUTINE BASIS_FINALISE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASIS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ALLOCATED(BASIS%INTERPOLATION_XI)) DEALLOCATE(BASIS%INTERPOLATION_XI)
      IF(ALLOCATED(BASIS%INTERPOLATION_TYPE)) DEALLOCATE(BASIS%INTERPOLATION_TYPE)
      IF(ALLOCATED(BASIS%INTERPOLATION_ORDER)) DEALLOCATE(BASIS%INTERPOLATION_ORDER)
      IF(ALLOCATED(BASIS%COLLAPSED_XI)) DEALLOCATE(BASIS%COLLAPSED_XI)
      IF(ALLOCATED(BASIS%NODE_AT_COLLAPSE)) DEALLOCATE(BASIS%NODE_AT_COLLAPSE)
      CALL BASIS_QUADRATURE_FINALISE(BASIS,ERR,ERROR,*999)
      IF(ALLOCATED(BASIS%NUMBER_OF_NODES_XIC)) DEALLOCATE(BASIS%NUMBER_OF_NODES_XIC)
      IF(ALLOCATED(BASIS%NUMBER_OF_DERIVATIVES)) DEALLOCATE(BASIS%NUMBER_OF_DERIVATIVES)
      IF(ALLOCATED(BASIS%NODE_POSITION_INDEX)) DEALLOCATE(BASIS%NODE_POSITION_INDEX)
      IF(ALLOCATED(BASIS%NODE_POSITION_INDEX_INV)) DEALLOCATE(BASIS%NODE_POSITION_INDEX_INV)
      IF(ALLOCATED(BASIS%DERIVATIVE_ORDER_INDEX)) DEALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX)
      IF(ALLOCATED(BASIS%DERIVATIVE_ORDER_INDEX_INV)) DEALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX_INV)
      IF(ALLOCATED(BASIS%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(BASIS%PARTIAL_DERIVATIVE_INDEX)
      IF(ALLOCATED(BASIS%ELEMENT_PARAMETER_INDEX)) DEALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX)
      IF(ALLOCATED(BASIS%ELEMENT_PARAMETER_INDEX_INV)) DEALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX_INV)
      IF(ALLOCATED(BASIS%LOCAL_LINE_XI_DIRECTION)) DEALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION)
      IF(ALLOCATED(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE)) DEALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE)
      IF(ALLOCATED(BASIS%NODE_NUMBERS_IN_LOCAL_LINE)) DEALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE)
      IF(ALLOCATED(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE)) DEALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE)
      IF(ALLOCATED(BASIS%LOCAL_FACE_XI_DIRECTION)) DEALLOCATE(BASIS%LOCAL_FACE_XI_DIRECTION)
      IF(ALLOCATED(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE)) DEALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE)
      IF(ALLOCATED(BASIS%NODE_NUMBERS_IN_LOCAL_FACE)) DEALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_FACE)
      IF(ALLOCATED(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE)) DEALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE)
      IF(ALLOCATED(BASIS%LOCAL_XI_NORMAL)) DEALLOCATE(BASIS%LOCAL_XI_NORMAL)
      IF(ASSOCIATED(BASIS%LINE_BASES)) DEALLOCATE(BASIS%LINE_BASES)
      IF(ASSOCIATED(BASIS%FACE_BASES)) DEALLOCATE(BASIS%FACE_BASES)
      IF(ASSOCIATED(BASIS%SUB_BASES)) DEALLOCATE(BASIS%SUB_BASES)      
      DEALLOCATE(BASIS)
    ENDIF
   
    CALL EXITS("BASIS_FINALISE")
    RETURN
999 CALL ERRORS("BASIS_FINALISE",ERR,ERROR)
    CALL EXITS("BASIS_FINALISE")
    RETURN 1
    
  END SUBROUTINE BASIS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a basis.
  SUBROUTINE BASIS_INITIALISE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASIS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      BASIS%USER_NUMBER=0
      BASIS%GLOBAL_NUMBER=0
      BASIS%FAMILY_NUMBER=0
      BASIS%BASIS_FINISHED=.FALSE.
      BASIS%HERMITE=.FALSE.
      BASIS%TYPE=0
      BASIS%NUMBER_OF_XI=0
      BASIS%NUMBER_OF_XI_COORDINATES=0
      BASIS%DEGENERATE=.FALSE.
      BASIS%NUMBER_OF_COLLAPSED_XI=0
      BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=0
      BASIS%NUMBER_OF_NODES=0
      BASIS%NUMBER_OF_ELEMENT_PARAMETERS=0
      BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES=0
      BASIS%NUMBER_OF_LOCAL_LINES=0
      BASIS%NUMBER_OF_LOCAL_FACES=0
      NULLIFY(BASIS%LINE_BASES)
      NULLIFY(BASIS%FACE_BASES)
      BASIS%NUMBER_OF_SUB_BASES=0
      NULLIFY(BASIS%SUB_BASES)
      NULLIFY(BASIS%PARENT_BASIS)
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
   
    CALL EXITS("BASIS_INITIALISE")
    RETURN
999 CALL ERRORS("BASIS_INITIALISE",ERR,ERROR)
    CALL EXITS("BASIS_INITIALISE")
    RETURN 1
    
  END SUBROUTINE BASIS_INITIALISE

  !
  !================================================================================================================================
  ! 

  !>Interpolates the appropriate partial derivative index of the element parameters at a gauss point for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !!>coordinate system with COORDINATE_INTERPOLATE_ADJUST. 
  FUNCTION BASIS_INTERPOLATE_GAUSS_DP(BASIS,PARTIAL_DERIV_INDEX,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER,ELEMENT_PARAMETERS,ERR,ERROR)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), INTENT(IN) :: QUADRATURE_SCHEME !<The quadrature scheme to use \see BASIS_ROUTINE_QuadratureSchemes
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The Gauss point number in the scheme to interpolte
    REAL(DP), INTENT(IN) :: ELEMENT_PARAMETERS(:) !<The element parameters to interpolate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_INTERPOLATE_GAUSS_DP
    !Local Variables
    INTEGER(INTG) :: ns
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: BASIS_QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BASIS_INTERPOLATE_GAUSS_DP",ERR,ERROR,*999)
    
    BASIS_INTERPOLATE_GAUSS_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(QUADRATURE_SCHEME>0.AND.QUADRATURE_SCHEME<=BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES) THEN
        BASIS_QUADRATURE_SCHEME=>BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(QUADRATURE_SCHEME)%PTR
        IF(ASSOCIATED(BASIS_QUADRATURE_SCHEME)) THEN
          IF(GAUSS_POINT_NUMBER>0.AND.GAUSS_POINT_NUMBER<=BASIS_QUADRATURE_SCHEME%NUMBER_OF_GAUSS) THEN
            IF(PARTIAL_DERIV_INDEX>0.AND.PARTIAL_DERIV_INDEX<=BASIS%NUMBER_OF_PARTIAL_DERIVATIVES) THEN
              DO ns=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                BASIS_INTERPOLATE_GAUSS_DP=BASIS_INTERPOLATE_GAUSS_DP+ &
                  & BASIS_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIV_INDEX,GAUSS_POINT_NUMBER)* &
                  & ELEMENT_PARAMETERS(ns)
              ENDDO !ns
            ELSE
              LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIV_INDEX,"*",ERR,ERROR))// &
                & " is invalid. It must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_PARTIAL_DERIVATIVES,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("The quadrature scheme has not been created",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The quadrature scheme type number of "//TRIM(NUMBER_TO_VSTRING(QUADRATURE_SCHEME,"*",ERR,ERROR))// &
          & " is invalid. It must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_INTERPOLATE_GAUSS_DP")
    RETURN
999 CALL ERRORS("BASIS_INTERPOLATE_GAUSS_DP",ERR,ERROR)
    CALL EXITS("BASIS_INTERPOLATE_GAUSS_DP")
  END FUNCTION BASIS_INTERPOLATE_GAUSS_DP

  !
  !================================================================================================================================
  !

  !>Interpolates the appropriate partial derivative index of the element local face parameters at a face gauss point for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !!>coordinate system with COORDINATE_INTERPOLATE_ADJUST. 
  FUNCTION BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP(BASIS,PARTIAL_DERIV_INDEX,QUADRATURE_SCHEME, &
    & LOCAL_FACE_NUMBER,GAUSS_POINT_NUMBER,FACE_PARAMETERS,ERR,ERROR)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), INTENT(IN) :: QUADRATURE_SCHEME !<The quadrature scheme to use \see BASIS_ROUTINE_QuadratureSchemes
    INTEGER(INTG), INTENT(IN) :: LOCAL_FACE_NUMBER !<The index number of the face to interpolate on
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The face Gauss point number in the scheme to interpolate
    REAL(DP), INTENT(IN) :: FACE_PARAMETERS(:) !<The face parameters to interpolate (in 3D coordinates)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP
    !Local Variables
    INTEGER(INTG) :: ns
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: BASIS_QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP",ERR,ERROR,*999)
    
    BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(QUADRATURE_SCHEME>0.AND.QUADRATURE_SCHEME<=BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES) THEN
        BASIS_QUADRATURE_SCHEME=>BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(QUADRATURE_SCHEME)%PTR
        IF(ASSOCIATED(BASIS_QUADRATURE_SCHEME)) THEN
          IF(BASIS%QUADRATURE%EVALUATE_FACE_GAUSS) THEN !Alternartively, can check whether scheme's face arrays are allocated?
            IF(LOCAL_FACE_NUMBER>0.AND.LOCAL_FACE_NUMBER<=BASIS%NUMBER_OF_LOCAL_FACES) THEN
              IF(GAUSS_POINT_NUMBER>0.AND.GAUSS_POINT_NUMBER<=BASIS_QUADRATURE_SCHEME%NUMBER_OF_FACE_GAUSS(LOCAL_FACE_NUMBER)) THEN
                IF(PARTIAL_DERIV_INDEX>0.AND.PARTIAL_DERIV_INDEX<=BASIS%NUMBER_OF_PARTIAL_DERIVATIVES) THEN
                  DO ns=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP=BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP+ &
                      & BASIS_QUADRATURE_SCHEME%FACE_GAUSS_BASIS_FNS(ns,PARTIAL_DERIV_INDEX,GAUSS_POINT_NUMBER,LOCAL_FACE_NUMBER)* &
                      & FACE_PARAMETERS(ns)
                  ENDDO !ns
                ELSE
                  LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIV_INDEX,"*",ERR,ERROR))// &
                    & " is invalid. It must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_PARTIAL_DERIVATIVES,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("The local face number index is invalid.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The face gauss interpolation scheme has not been created",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The quadrature scheme has not been created",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The quadrature scheme type number of "//TRIM(NUMBER_TO_VSTRING(QUADRATURE_SCHEME,"*",ERR,ERROR))// &
          & " is invalid. It must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP")
    RETURN
999 CALL ERRORS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP",ERR,ERROR)
    CALL EXITS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP")
  END FUNCTION BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP

  !
  !================================================================================================================================
  !

  !>Interpolates the appropriate partial derivative index of the element parameters at position XI for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !>coordinate system with COORDINATE_INTERPOLATE_ADJUST. Note for simplex basis functions the XI coordinates should
  !>exclude the last area coordinate.
  FUNCTION BASIS_INTERPOLATE_XI_DP(BASIS,PARTIAL_DERIV_INDEX,XI,ELEMENT_PARAMETERS,ERR,ERROR)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: XI(:) !<The Xi position to interpolate the basis function at
    REAL(DP), INTENT(IN) :: ELEMENT_PARAMETERS(:) !<The element parameters to interpolate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_INTERPOLATE_XI_DP
    !Local Variables
    INTEGER(INTG) :: nn,nk,ns
    REAL(DP) :: XIL(SIZE(XI,1)+1)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BASIS_INTERPOLATE_XI_DP",ERR,ERROR,*999)
    
    BASIS_INTERPOLATE_XI_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      SELECT CASE(BASIS%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        ns=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
            ns=ns+1
            BASIS_INTERPOLATE_XI_DP=BASIS_INTERPOLATE_XI_DP+ &
              & BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,PARTIAL_DERIV_INDEX,XI,ERR,ERROR)* &
              & ELEMENT_PARAMETERS(ns)
          ENDDO !nk
        ENDDO !nn
        IF(ERR/=0) GOTO 999
      CASE(BASIS_SIMPLEX_TYPE)
        !Create the area coordinates from the xi coordinates
        XIL(1:SIZE(XI,1))=1.0_DP-XI
        XIL(SIZE(XI,1)+1)=SUM(XI)-(SIZE(XI,1)-1.0_DP)
        ns=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          ns=ns+1
          BASIS_INTERPOLATE_XI_DP=BASIS_INTERPOLATE_XI_DP+ &
            & BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,PARTIAL_DERIV_INDEX,XIL,ERR,ERROR)* &
            & ELEMENT_PARAMETERS(ns)
        ENDDO !nn
        IF(ERR/=0) GOTO 999
      CASE DEFAULT
        LOCAL_ERROR="Basis type "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))//" is invalid or not implemented"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_INTERPOLATE_XI_DP")
    RETURN
999 CALL ERRORS("BASIS_INTERPOLATE_XI_DP",ERR,ERROR)
    CALL EXITS("BASIS_INTERPOLATE_XI_DP")
  END FUNCTION BASIS_INTERPOLATE_XI_DP
  
  !
  !================================================================================================================================
  !
  
  !>Gets/changes the interpolation type in each xi directions for a basis identified by a pointer.
  SUBROUTINE BASIS_INTERPOLATION_XI_GET(BASIS,INTERPOLATION_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to get the interpolation xi
    INTEGER(INTG), INTENT(OUT) :: INTERPOLATION_XI(:) !<On return, the interpolation xi parameters for each Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_INTERPOLATION_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(SIZE(INTERPOLATION_XI,1)>=SIZE(BASIS%INTERPOLATION_XI,1)) THEN
          INTERPOLATION_XI=BASIS%INTERPOLATION_XI
        ELSE
          LOCAL_ERROR="The size of INTERPOLATION_XI is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(INTERPOLATION_XI,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(BASIS%INTERPOLATION_XI,1),"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Basis has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_INTERPOLATION_XI_GET")
    RETURN
999 CALL ERRORS("BASIS_INTERPOLATION_XI_GET",ERR,ERROR)
    CALL EXITS("BASIS_INTERPOLATION_XI_GET")
    RETURN
  END SUBROUTINE BASIS_INTERPOLATION_XI_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type in each xi directions where the basis is identified by user number. \see OPENCMISS::CMISSBasisInterpolationXiSet
  SUBROUTINE BASIS_INTERPOLATION_XI_SET_NUMBER(USER_NUMBER,INTERPOLATION_XI,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to set the interpolation xi
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_XI(:) !<The interpolation xi parameters for each Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_INTERPOLATION_XI_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_SET(BASIS,INTERPOLATION_XI,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_INTERPOLATION_XI_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_INTERPOLATION_XI_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_INTERPOLATION_XI_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_INTERPOLATION_XI_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type in each xi directions for a basis identified by a pointer. \see OPENCMISS::CMISSBasisInterpolationXiSet
  SUBROUTINE BASIS_INTERPOLATION_XI_SET_PTR(BASIS,INTERPOLATION_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set the interpolation xi
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_XI(:) !<The interpolation xi parameters for each Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni,LAST_INTERP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_INTERPOLATION_XI_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        IF(SIZE(INTERPOLATION_XI,1)==BASIS%NUMBER_OF_XI) THEN
          !Check the input values
          SELECT CASE(BASIS%TYPE)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            DO ni=1,BASIS%NUMBER_OF_XI
              SELECT CASE(INTERPOLATION_XI(ni))
              CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
                & BASIS_CUBIC_HERMITE_INTERPOLATION,BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
                !Do nothing
              CASE DEFAULT
                LOCAL_ERROR="Interpolation xi value "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_XI(ni),"*",ERR,ERROR))// &
                  & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid for a Lagrange-Hermite TP basis"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !ni
          CASE(BASIS_SIMPLEX_TYPE)
            LAST_INTERP=INTERPOLATION_XI(1)
            DO ni=1,BASIS%NUMBER_OF_XI
              SELECT CASE(INTERPOLATION_XI(ni))
              CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,BASIS_CUBIC_SIMPLEX_INTERPOLATION)
                IF(INTERPOLATION_XI(ni)/=LAST_INTERP) THEN
                  CALL FLAG_ERROR("The interpolation xi value must be the same for all xi directions for a simplex basis", &
                    & ERR,ERROR,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="Interpolation xi value "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_XI(ni),"*",ERR,ERROR))// &
                  & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid for a simplex basis"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !ni
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid basis type or not implemented",ERR,ERROR,*999)
          END SELECT
          !Set the interpolation xi
          BASIS%INTERPOLATION_XI=INTERPOLATION_XI
        ELSE
          LOCAL_ERROR="The size of the interpolation xi array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(INTERPOLATION_XI,1),"*",ERR,ERROR))//") does not match the number of xi directions ("// &
            & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//") for basis number "// &
            & TRIM(NUMBER_TO_VSTRING(BASIS%USER_NUMBER,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_INTERPOLATION_XI_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_INTERPOLATION_XI_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_INTERPOLATION_XI_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_INTERPOLATION_XI_SET_PTR

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Lagrange-Hermite tensor product basis that has already been allocated BASIS_CREATE_START
  !> \see BASIS_ROUTINES::BASIS_CREATE_START
  SUBROUTINE BASIS_LHTP_BASIS_CREATE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MAX_NUM_NODES,NUMBER_OF_DERIVATIVES,ni,ni1,ni2,ni3,nk,nn,nnl,ns,OLD_NUMBER_OF_DERIVATIVES, &
      & POSITION(4),MAXIMUM_NODE_EXTENT(3),COLLAPSED_XI(3),NUMBER_OF_NODES,NUMBER_OF_LOCAL_LINES,NODE_COUNT, &
      & SPECIAL_NODE_COUNT,NODES_IN_LINE(4),nnf,NUMBER_OF_LOCAL_FACES,LOCAL_NODE_COUNT,ef, NUMBER_OF_ELEMENT_FACES
    INTEGER, TARGET :: nn1,nn2,nn3,nn4
    LOGICAL, ALLOCATABLE :: NODE_AT_COLLAPSE(:)
    !move to types.f90
    TYPE INTG_POINTER
       INTEGER, POINTER :: a1, a2, a3, a4
    END TYPE INTG_POINTER
    TYPE(INTG_POINTER) :: ARGLIST
    INTEGER, POINTER :: np1, np2, np3, np4
    LOGICAL :: AT_COLLAPSE,FIRST_COLLAPSED_POSITION,PROCESS_NODE
    
    CALL ENTERS("BASIS_LHTP_BASIS_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
        BASIS%NUMBER_OF_XI_COORDINATES=BASIS%NUMBER_OF_XI
        BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=BASIS%NUMBER_OF_XI_COORDINATES**2+2
        ALLOCATE(BASIS%INTERPOLATION_TYPE(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate the interpolation type array",ERR,ERROR,*999)
        ALLOCATE(BASIS%INTERPOLATION_ORDER(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation order array",ERR,ERROR,*999)
        ALLOCATE(BASIS%NUMBER_OF_NODES_XIC(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes xic array",ERR,ERROR,*999)
        NUMBER_OF_NODES=1
        MAX_NUM_NODES=0
        BASIS%DEGENERATE=.FALSE.
        BASIS%NUMBER_OF_COLLAPSED_XI=0
        DO ni=1,BASIS%NUMBER_OF_XI
          SELECT CASE(BASIS%INTERPOLATION_XI(ni))
          CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(ni)=BASIS_LAGRANGE_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(ni)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(ni)=2            
          CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(ni)=BASIS_LAGRANGE_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(ni)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(ni)=3
          CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(ni)=BASIS_LAGRANGE_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(ni)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(ni)=4
          CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(ni)=BASIS_HERMITE_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(ni)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(ni)=2
          CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(ni)=BASIS_HERMITE_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(ni)=BASIS_QUADRATIC1_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(ni)=2
          CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(ni)=BASIS_HERMITE_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(ni)=BASIS_QUADRATIC2_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(ni)=2
          CASE DEFAULT 
            CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
          END SELECT
          IF(BASIS%COLLAPSED_XI(ni)==BASIS_XI_COLLAPSED) THEN
            BASIS%NUMBER_OF_COLLAPSED_XI=BASIS%NUMBER_OF_COLLAPSED_XI+1
            COLLAPSED_XI(BASIS%NUMBER_OF_COLLAPSED_XI)=ni
            BASIS%DEGENERATE=.TRUE.
          ENDIF
          NUMBER_OF_NODES=NUMBER_OF_NODES*BASIS%NUMBER_OF_NODES_XIC(ni)
          IF(BASIS%NUMBER_OF_NODES_XIC(ni)>MAX_NUM_NODES) MAX_NUM_NODES=BASIS%NUMBER_OF_NODES_XIC(ni)
        ENDDO !ni
        !If a degenerate (collapsed) basis recalculate the number of nodes from the maximum posible number of nodes
        IF(BASIS%DEGENERATE) THEN
          ALLOCATE(NODE_AT_COLLAPSE(NUMBER_OF_NODES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate at collapse",ERR,ERROR,*999)
          POSITION=1
          BASIS%NUMBER_OF_NODES=0
          !Loop over the maximum number of nodes which is currently set for the basis
          DO nn=1,NUMBER_OF_NODES
            AT_COLLAPSE=.FALSE.
            DO ni=1,BASIS%NUMBER_OF_XI
              IF(BASIS%COLLAPSED_XI(ni)==BASIS_COLLAPSED_AT_XI0.AND.POSITION(ni)==1.OR. &
                & BASIS%COLLAPSED_XI(ni)==BASIS_COLLAPSED_AT_XI1.AND.POSITION(ni)==BASIS%NUMBER_OF_NODES_XIC(ni)) &
                & AT_COLLAPSE=.TRUE.
            ENDDO !ni
            IF(AT_COLLAPSE) THEN
              IF(BASIS%NUMBER_OF_COLLAPSED_XI==1) THEN
                FIRST_COLLAPSED_POSITION=POSITION(COLLAPSED_XI(1))==1
              ELSE
                FIRST_COLLAPSED_POSITION=(POSITION(COLLAPSED_XI(1))==1).AND.(POSITION(COLLAPSED_XI(2))==1)
              ENDIF
              IF(FIRST_COLLAPSED_POSITION) THEN
                BASIS%NUMBER_OF_NODES=BASIS%NUMBER_OF_NODES+1
                NODE_AT_COLLAPSE(BASIS%NUMBER_OF_NODES)=.TRUE.
              ENDIF
            ELSE
              BASIS%NUMBER_OF_NODES=BASIS%NUMBER_OF_NODES+1
              NODE_AT_COLLAPSE(BASIS%NUMBER_OF_NODES)=.FALSE.
            ENDIF
            POSITION(1)=POSITION(1)+1
            DO ni=1,BASIS%NUMBER_OF_XI
              IF(POSITION(ni)>BASIS%NUMBER_OF_NODES_XIC(ni)) THEN
                POSITION(ni)=1
                POSITION(ni+1)=POSITION(ni+1)+1
              ENDIF
            ENDDO !ni
          ENDDO !nn
          ALLOCATE(BASIS%NODE_AT_COLLAPSE(BASIS%NUMBER_OF_NODES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis node at collapse",ERR,ERROR,*999)
          BASIS%NODE_AT_COLLAPSE(1:BASIS%NUMBER_OF_NODES)=NODE_AT_COLLAPSE(1:BASIS%NUMBER_OF_NODES)
          DEALLOCATE(NODE_AT_COLLAPSE)
        ELSE        
          BASIS%NUMBER_OF_NODES=NUMBER_OF_NODES
          ALLOCATE(BASIS%NODE_AT_COLLAPSE(BASIS%NUMBER_OF_NODES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis node at collapse",ERR,ERROR,*999)
          BASIS%NODE_AT_COLLAPSE=.FALSE.
          COLLAPSED_XI(1)=1
        ENDIF
        
        ALLOCATE(BASIS%NODE_POSITION_INDEX(BASIS%NUMBER_OF_NODES,BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODE_POSITION_INDEX",ERR,ERROR,*999)
        SELECT CASE(BASIS%NUMBER_OF_XI_COORDINATES)
        CASE(1)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,1,1,1),STAT=ERR)
        CASE(2)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,1,1),STAT=ERR)
        CASE(3)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES,1),STAT=ERR)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
        END SELECT
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODE_POSITION_INDEX_INV",ERR,ERROR,*999)
        BASIS%NODE_POSITION_INDEX_INV=0
        
        !Determine the node position index and it's inverse
        POSITION=1
        nn=0
        FIRST_COLLAPSED_POSITION=.TRUE.
        DO nn1=1,NUMBER_OF_NODES
          AT_COLLAPSE=.FALSE.
          IF(BASIS%DEGENERATE) THEN
            DO ni=1,BASIS%NUMBER_OF_XI
              IF(BASIS%COLLAPSED_XI(ni)==BASIS_COLLAPSED_AT_XI0.AND.POSITION(ni)==1.OR. &
                & BASIS%COLLAPSED_XI(ni)==BASIS_COLLAPSED_AT_XI1.AND.POSITION(ni)==BASIS%NUMBER_OF_NODES_XIC(ni)) &
                & AT_COLLAPSE=.TRUE.
            ENDDO !ni
            IF(BASIS%NUMBER_OF_COLLAPSED_XI==1) THEN
              FIRST_COLLAPSED_POSITION=POSITION(COLLAPSED_XI(1))==1
            ELSE
              FIRST_COLLAPSED_POSITION=(POSITION(COLLAPSED_XI(1))==1).AND.(POSITION(COLLAPSED_XI(2))==1)
            ENDIF
          ENDIF
          PROCESS_NODE=(AT_COLLAPSE.AND.FIRST_COLLAPSED_POSITION).OR.(.NOT.AT_COLLAPSE)
          IF(PROCESS_NODE) THEN
            nn=nn+1
            BASIS%NODE_POSITION_INDEX(nn,:)=POSITION(1:BASIS%NUMBER_OF_XI)
!!Should the inverse of the node position index be adjusted so that the collapsed positions point to the collapsed nn?
!!At the moment this is not the case and they are just set to zero.
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%NODE_POSITION_INDEX_INV(BASIS%NODE_POSITION_INDEX(nn,1),1,1,1)=nn
            CASE(2)
              BASIS%NODE_POSITION_INDEX_INV(BASIS%NODE_POSITION_INDEX(nn,1),BASIS%NODE_POSITION_INDEX(nn,2),1,1)=nn
            CASE(3)
              BASIS%NODE_POSITION_INDEX_INV(BASIS%NODE_POSITION_INDEX(nn,1),BASIS%NODE_POSITION_INDEX(nn,2), &
                & BASIS%NODE_POSITION_INDEX(nn,3),1)=nn
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of Xi directions",ERR,ERROR,*999)
            END SELECT
          ENDIF
          POSITION(1)=POSITION(1)+1
          DO ni=1,BASIS%NUMBER_OF_XI
            IF(POSITION(ni)>BASIS%NUMBER_OF_NODES_XIC(ni)) THEN
              POSITION(ni)=1
              POSITION(ni+1)=POSITION(ni+1)+1
            ENDIF
          ENDDO !ni
        ENDDO !nn1
        !Calculate the maximum number of derivatives and number of element
        !parameters
        BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES=-1
        BASIS%NUMBER_OF_ELEMENT_PARAMETERS=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          NUMBER_OF_DERIVATIVES=1
          DO ni=1,BASIS%NUMBER_OF_XI
            IF((.NOT.BASIS%NODE_AT_COLLAPSE(nn).OR.BASIS%COLLAPSED_XI(ni)==BASIS_NOT_COLLAPSED).AND. &
              & BASIS%INTERPOLATION_TYPE(ni)==BASIS_HERMITE_INTERPOLATION.AND. &
              & (BASIS%INTERPOLATION_ORDER(ni)==BASIS_CUBIC_INTERPOLATION_ORDER.OR. &
              & (BASIS%NODE_POSITION_INDEX(nn,ni)==1.AND. &
              & BASIS%INTERPOLATION_ORDER(ni)==BASIS_QUADRATIC2_INTERPOLATION_ORDER).OR. &
              & (BASIS%NODE_POSITION_INDEX(nn,ni)==2.AND. &
              & BASIS%INTERPOLATION_ORDER(ni)==BASIS_QUADRATIC1_INTERPOLATION_ORDER))) THEN
              !Derivative in this direction
              NUMBER_OF_DERIVATIVES=NUMBER_OF_DERIVATIVES*2
            ENDIF
          ENDDO !ni
          BASIS%NUMBER_OF_ELEMENT_PARAMETERS=BASIS%NUMBER_OF_ELEMENT_PARAMETERS+NUMBER_OF_DERIVATIVES
          IF(NUMBER_OF_DERIVATIVES>BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES) BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES= &
            & NUMBER_OF_DERIVATIVES
        ENDDO !nn
        !Now set up the number of derivatives and derivative order index
        ALLOCATE(BASIS%NUMBER_OF_DERIVATIVES(BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NUMBER_OF_DERIVATIVES",ERR,ERROR,*999)
        ALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES, &
          & BASIS%NUMBER_OF_XI),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DERIVATIVE_ORDER_INDEX",ERR,ERROR,*999)
        ALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX_INV(FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV, &
          & BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DERIVATIVE_ORDER_INDEX_INV",ERR,ERROR,*999)
        ALLOCATE(BASIS%PARTIAL_DERIVATIVE_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate PARTIAL_DERIVATIVE_INDEX",ERR,ERROR,*999)
        ALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT_PARAMETER_INDEX",ERR,ERROR,*999)
        ALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX_INV(2,BASIS%NUMBER_OF_ELEMENT_PARAMETERS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT_PARAMETER_INDEX_INV",ERR,ERROR,*999)
        !Set the derivative order index and its inverse, the element parameter index and the partial derivative index.
        ns=0
        BASIS%DERIVATIVE_ORDER_INDEX=0
        BASIS%DERIVATIVE_ORDER_INDEX_INV=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          BASIS%NUMBER_OF_DERIVATIVES(nn)=1
          DO ni=1,BASIS%NUMBER_OF_XI
            IF((.NOT.BASIS%NODE_AT_COLLAPSE(nn).OR.BASIS%COLLAPSED_XI(ni)==BASIS_NOT_COLLAPSED).AND. &
              & BASIS%INTERPOLATION_TYPE(ni)==BASIS_HERMITE_INTERPOLATION.AND. &
              & (BASIS%INTERPOLATION_ORDER(ni)==BASIS_CUBIC_INTERPOLATION_ORDER.OR. &
              & (BASIS%NODE_POSITION_INDEX(nn,ni)==1.AND. &
              & BASIS%INTERPOLATION_ORDER(ni)==BASIS_QUADRATIC2_INTERPOLATION_ORDER).OR. &
              & (BASIS%NODE_POSITION_INDEX(nn,ni)==2.AND. &
              & BASIS%INTERPOLATION_ORDER(ni)==BASIS_QUADRATIC1_INTERPOLATION_ORDER))) THEN
              OLD_NUMBER_OF_DERIVATIVES=BASIS%NUMBER_OF_DERIVATIVES(nn)
              BASIS%NUMBER_OF_DERIVATIVES(nn)=BASIS%NUMBER_OF_DERIVATIVES(nn)*2
              DO nk=1,OLD_NUMBER_OF_DERIVATIVES
                BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,ni)=NO_PART_DERIV
                BASIS%DERIVATIVE_ORDER_INDEX(OLD_NUMBER_OF_DERIVATIVES+nk,nn,ni)=FIRST_PART_DERIV
                DO ni2=1,ni-1
                  BASIS%DERIVATIVE_ORDER_INDEX(OLD_NUMBER_OF_DERIVATIVES+nk,nn,ni2)=BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,ni2)
                ENDDO !ni
              ENDDO !nk
            ELSE
              DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,ni)=1
              ENDDO !nk
            ENDIF
          ENDDO !ni

          DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
            ns=ns+1
            BASIS%ELEMENT_PARAMETER_INDEX(nk,nn)=ns
            BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ns)=nn
            BASIS%ELEMENT_PARAMETER_INDEX_INV(2,ns)=nk
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%DERIVATIVE_ORDER_INDEX_INV(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,1),1,1,nn)=nk
              SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,1))
              CASE(NO_PART_DERIV)
                BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=NO_PART_DERIV
              CASE(FIRST_PART_DERIV)
                BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1
              CASE DEFAULT
                CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
              END SELECT
            CASE(2)
              BASIS%DERIVATIVE_ORDER_INDEX_INV(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,1), &
                & BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,2),1,nn)=nk            
              SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,1))
              CASE(NO_PART_DERIV)
                SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,2))
                CASE(NO_PART_DERIV)
                  BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=NO_PART_DERIV
                CASE(FIRST_PART_DERIV)
                  BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S2
                CASE DEFAULT
                  CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                END SELECT
              CASE(FIRST_PART_DERIV)
                SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,2))
                CASE(NO_PART_DERIV)
                  BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1
                CASE(FIRST_PART_DERIV)
                  BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1_S2
                CASE DEFAULT
                  CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                END SELECT
              CASE DEFAULT
                CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
              END SELECT
            CASE(3)
              BASIS%DERIVATIVE_ORDER_INDEX_INV(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,1), &
                & BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,2),BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,3),nn)=nk           
              SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,1))
              CASE(NO_PART_DERIV)
                SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,2))
                CASE(NO_PART_DERIV)
                  SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,3))
                  CASE(NO_PART_DERIV)                
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=NO_PART_DERIV
                  CASE(FIRST_PART_DERIV)
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S3
                  CASE DEFAULT
                    CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                  END SELECT
                CASE(FIRST_PART_DERIV)
                  SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,3))
                  CASE(NO_PART_DERIV)                
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S2
                  CASE(FIRST_PART_DERIV)
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S2_S3
                  CASE DEFAULT
                    CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                  END SELECT
                CASE DEFAULT
                  CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                END SELECT
              CASE(FIRST_PART_DERIV)
                SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,2))
                CASE(NO_PART_DERIV)
                  SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,3))
                  CASE(NO_PART_DERIV)                
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1
                  CASE(FIRST_PART_DERIV)
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1_S3
                  CASE DEFAULT
                    CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                  END SELECT
                CASE(FIRST_PART_DERIV)
                  SELECT CASE(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn,3))
                  CASE(NO_PART_DERIV)                
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1_S2
                  CASE(FIRST_PART_DERIV)
                    BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)=PART_DERIV_S1_S2_S3
                  CASE DEFAULT
                    CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                  END SELECT
                CASE DEFAULT
                  CALL FLAG_ERROR("Invalid derivative order index",ERR,ERROR,*999)
                END SELECT
              CASE DEFAULT
                CALL FLAG_ERROR("Invalid derivative order index.",ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of Xi direcions.",ERR,ERROR,*999)
            END SELECT
          ENDDO !nk
        ENDDO !nn

        !Set up the line information
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          NUMBER_OF_LOCAL_LINES=1
          BASIS%NUMBER_OF_LOCAL_LINES=1
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local line.",ERR,ERROR,*999)
          BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=BASIS%NUMBER_OF_NODES_XIC(1)
          ALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line xi direction.",ERR,ERROR,*999)
          BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(BASIS%NUMBER_OF_NODES_XIC(1),NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local line.",ERR,ERROR,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(BASIS%NUMBER_OF_NODES_XIC(1),NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local line.",ERR,ERROR,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          DO nn2=1,BASIS%NUMBER_OF_NODES_XIC(1)
            DO nn1=1,BASIS%NUMBER_OF_NODES
              IF(BASIS%NODE_POSITION_INDEX(nn1,1)==nn2) THEN
                BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nn2,1)=nn1
                DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn2)
                  IF(BASIS%DERIVATIVE_ORDER_INDEX(nk,nn2,1)==FIRST_PART_DERIV) THEN
                    BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nn2,1)=nk
                    EXIT
                  ENDIF
                ENDDO !nk                
                EXIT
              ENDIF
            ENDDO !nn2
          ENDDO !nn1
        CASE(2)
          !Determine the maximum node extent of the basis
          MAXIMUM_NODE_EXTENT(1)=MAXVAL(BASIS%NODE_POSITION_INDEX(:,1))
          MAXIMUM_NODE_EXTENT(2)=MAXVAL(BASIS%NODE_POSITION_INDEX(:,2))
          !Allocate and calculate the lines
          NUMBER_OF_LOCAL_LINES=4-BASIS%NUMBER_OF_COLLAPSED_XI
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local line.",ERR,ERROR,*999)
          BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE=0
          ALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line xi direction.",ERR,ERROR,*999)
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local line",ERR,ERROR,*999)
          BASIS%NODE_NUMBERS_IN_LOCAL_LINE=0
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local line.",ERR,ERROR,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(BASIS%LOCAL_XI_NORMAL(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local xi normal.",ERR,ERROR,*999)
          !Find the lines
          BASIS%NUMBER_OF_LOCAL_LINES=0
          DO ni1=1,2
            ni2=OTHER_XI_DIRECTIONS2(ni1)
            !We are looking for lines in the ni1 direction from the direction of ni1=0
            !Loop over the element extremes in the ni2 direction
            DO nn2=1,MAXIMUM_NODE_EXTENT(ni2),MAXIMUM_NODE_EXTENT(ni2)-1
              NODE_COUNT=0
              SPECIAL_NODE_COUNT=0
              NODES_IN_LINE=0
              DO nn1=1,BASIS%NUMBER_OF_NODES
                IF(BASIS%COLLAPSED_XI(ni1)/=BASIS_NOT_COLLAPSED) THEN
                  !The current xi direction, ni1, is in a degenerate plane
                  IF(BASIS%COLLAPSED_XI(ni2)==BASIS_XI_COLLAPSED) THEN
                    !The other xi direction is collapsed (must be the case)
                    IF(BASIS%COLLAPSED_XI(ni1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                      IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.OR.BASIS%NODE_POSITION_INDEX(nn1,ni1)==1) THEN
                        NODE_COUNT=NODE_COUNT+1
                        NODES_IN_LINE(NODE_COUNT)=nn1
                      ENDIF
                    ELSE !Collapsed at the xi=1 end
                      IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2) THEN
                        NODE_COUNT=NODE_COUNT+1
                        NODES_IN_LINE(NODE_COUNT)=nn1
                      ELSE IF(BASIS%NODE_POSITION_INDEX(nn1,ni1)==MAXIMUM_NODE_EXTENT(ni1)) THEN
                        IF(ni1<2) THEN !Special case - put the collapsed node at the end of the line
                          SPECIAL_NODE_COUNT=SPECIAL_NODE_COUNT+1
                          NODES_IN_LINE(MAXIMUM_NODE_EXTENT(ni1))=nn1
                        ELSE
                          NODE_COUNT=NODE_COUNT+1
                          NODES_IN_LINE(NODE_COUNT)=nn1
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    !The current xi direction must be collapsed
                    IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2) THEN
                      NODE_COUNT=NODE_COUNT+1
                      NODES_IN_LINE(NODE_COUNT)=nn1
                    ENDIF
                  ENDIF
                ELSE
                  !The current xi direction, ni1, is not involved in any collapsed (degenerate) planes
                  IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2) THEN
                    NODE_COUNT=NODE_COUNT+1
                    NODES_IN_LINE(NODE_COUNT)=nn1
                  ENDIF
                ENDIF
              ENDDO !nn1
              IF((NODE_COUNT+SPECIAL_NODE_COUNT)>1) THEN !More than one node so it is a proper line 
                BASIS%NUMBER_OF_LOCAL_LINES=BASIS%NUMBER_OF_LOCAL_LINES+1
                BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES)=NODE_COUNT+SPECIAL_NODE_COUNT
                BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES), &
                  & BASIS%NUMBER_OF_LOCAL_LINES)=NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES))
                DO nnl=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES)
                  DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,BASIS%NUMBER_OF_LOCAL_LINES))
                    IF(BASIS%DERIVATIVE_ORDER_INDEX(nk,BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,BASIS%NUMBER_OF_LOCAL_LINES),ni1)== &
                      & FIRST_PART_DERIV) THEN
                      BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nnl,BASIS%NUMBER_OF_LOCAL_LINES)=nk
                      EXIT
                    ENDIF
                  ENDDO !nk
                ENDDO !nnl
                BASIS%LOCAL_LINE_XI_DIRECTION(BASIS%NUMBER_OF_LOCAL_LINES)=ni1
                IF(nn2==1) THEN
                  BASIS%LOCAL_XI_NORMAL(BASIS%NUMBER_OF_LOCAL_LINES)=-ni2
                ELSE
                  BASIS%LOCAL_XI_NORMAL(BASIS%NUMBER_OF_LOCAL_LINES)=ni2
                ENDIF
              ENDIF
            ENDDO !nn2
          ENDDO !ni1
        CASE(3)
          !Determine the maximum node extent of the basis
          MAXIMUM_NODE_EXTENT(1)=MAXVAL(BASIS%NODE_POSITION_INDEX(:,1))
          MAXIMUM_NODE_EXTENT(2)=MAXVAL(BASIS%NODE_POSITION_INDEX(:,2))
          MAXIMUM_NODE_EXTENT(3)=MAXVAL(BASIS%NODE_POSITION_INDEX(:,3))
          !Allocate and calculate the lines
          IF(BASIS%NUMBER_OF_COLLAPSED_XI==1) THEN
            NUMBER_OF_LOCAL_LINES=9
            NUMBER_OF_LOCAL_FACES=5
            BASIS%NUMBER_OF_LOCAL_FACES=5
          ELSE IF(BASIS%NUMBER_OF_COLLAPSED_XI==2) THEN
            NUMBER_OF_LOCAL_LINES=8
            NUMBER_OF_LOCAL_FACES=5
            BASIS%NUMBER_OF_LOCAL_FACES=5
          ELSE
            NUMBER_OF_LOCAL_LINES=12
            NUMBER_OF_LOCAL_FACES=6
            BASIS%NUMBER_OF_LOCAL_FACES=6
          ENDIF

          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local line",ERR,ERROR,*999)
          BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE=0

          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local face",ERR,ERROR,*999)
          BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE=0

          ALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION(NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line xi direction",ERR,ERROR,*999)

          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local line",ERR,ERROR,*999)
          BASIS%NODE_NUMBERS_IN_LOCAL_LINE=0

          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local line",ERR,ERROR,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV

          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC)**2,NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local face",ERR,ERROR,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE=NO_PART_DERIV

          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(MAX(MAXIMUM_NODE_EXTENT(2)*MAXIMUM_NODE_EXTENT(3), &
                               & MAXIMUM_NODE_EXTENT(3)*MAXIMUM_NODE_EXTENT(1), &
                               & MAXIMUM_NODE_EXTENT(2)*MAXIMUM_NODE_EXTENT(1)),NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local face",ERR,ERROR,*999)
          BASIS%NODE_NUMBERS_IN_LOCAL_FACE=0
          ALLOCATE(BASIS%LOCAL_XI_NORMAL(NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local xi normal.",ERR,ERROR,*999)
          
          ALLOCATE(BASIS%LOCAL_FACE_XI_DIRECTION(NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local face xi direction",ERR,ERROR,*999)
          
          !Find the lines and faces
          BASIS%NUMBER_OF_LOCAL_LINES=0
          DO ni1=1,3
            ni2=OTHER_XI_DIRECTIONS3(ni1,2,1)
            ni3=OTHER_XI_DIRECTIONS3(ni1,3,1)
            
           ! DO nn3=1,MAXIMUM_NODE_EXTENT(ni3)
            !  DO nn2=1,MAXIMUM_NODE_EXTENT(ni2)
                 
            !We are looking for lines going in the ni1 direction, starting from ni1=0.
            DO nn3=1,MAXIMUM_NODE_EXTENT(ni3),MAXIMUM_NODE_EXTENT(ni3)-1 
              DO nn2=1,MAXIMUM_NODE_EXTENT(ni2),MAXIMUM_NODE_EXTENT(ni2)-1
                NODE_COUNT=0
                SPECIAL_NODE_COUNT=0
                NODES_IN_LINE=0
                ! Iterate over nodes in the line of interest
                DO nn1=1,BASIS%NUMBER_OF_NODES
                  IF(BASIS%COLLAPSED_XI(ni1)/=BASIS_NOT_COLLAPSED) THEN
                    !The current xi direction, ni1, is involved in a collapsed (degenerate) plane 
                    IF(BASIS%COLLAPSED_XI(ni2)==BASIS_XI_COLLAPSED.AND.BASIS%COLLAPSED_XI(ni3)==BASIS_XI_COLLAPSED) THEN
                      !Both of the other two xi directions are collapsed
                      IF(BASIS%COLLAPSED_XI(ni1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                        IF((BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.OR.BASIS%NODE_POSITION_INDEX(nn1,ni1)==1).AND. &
                          & (BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3.OR.BASIS%NODE_POSITION_INDEX(nn1,ni1)==1)) THEN
                          NODE_COUNT=NODE_COUNT+1
                          NODES_IN_LINE(NODE_COUNT)=nn1
                        ENDIF
                      ELSE !Collapsed at the xi=1 end
                        IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                          NODE_COUNT=NODE_COUNT+1
                          NODES_IN_LINE(NODE_COUNT)=nn1
                        ELSE IF(BASIS%NODE_POSITION_INDEX(nn1,ni1)==MAXIMUM_NODE_EXTENT(ni1)) THEN
                          IF(ni1<3) THEN !Special case - put the collapsed node at the end of the line
                            SPECIAL_NODE_COUNT=SPECIAL_NODE_COUNT+1
                            NODES_IN_LINE(MAXIMUM_NODE_EXTENT(ni1))=nn1
                          ELSE
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE
                      IF(BASIS%COLLAPSED_XI(ni2)==BASIS_XI_COLLAPSED) THEN
                        !The other ni2 xi direction is collapsed
                        IF(BASIS%COLLAPSED_XI(ni1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                          IF((BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.OR.BASIS%NODE_POSITION_INDEX(nn1,ni1)==1).AND. &
                            & BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ENDIF
                        ELSE IF(BASIS%COLLAPSED_XI(ni1)==BASIS_COLLAPSED_AT_XI1) THEN !Collapsed at the xi=1 end
                          IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ELSE IF(BASIS%NODE_POSITION_INDEX(nn1,ni1)==MAXIMUM_NODE_EXTENT(ni1).AND. &
                            & BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                            IF(ni1<ni2) THEN !Special case - put the collapsed node at the end of the line
                              SPECIAL_NODE_COUNT=SPECIAL_NODE_COUNT+1
                              NODES_IN_LINE(MAXIMUM_NODE_EXTENT(ni1))=nn1
                            ELSE
                              NODE_COUNT=NODE_COUNT+1
                              NODES_IN_LINE(NODE_COUNT)=nn1
                            ENDIF
                          ENDIF
                        ELSE
                          !Not collapsed at a xi end
                          IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ENDIF
                        ENDIF
                      ELSE IF(BASIS%COLLAPSED_XI(ni3)==BASIS_XI_COLLAPSED) THEN
                        !The other ni3 xi direction is collapsed
                        IF(BASIS%COLLAPSED_XI(ni1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                          IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND. &
                            & (BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3.OR.BASIS%NODE_POSITION_INDEX(nn1,ni1)==1)) THEN
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ENDIF
                        ELSE IF(BASIS%COLLAPSED_XI(ni1)==BASIS_COLLAPSED_AT_XI1) THEN !Collapsed at the xi=1 end
                          IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ELSE IF(BASIS%NODE_POSITION_INDEX(nn1,ni1)==MAXIMUM_NODE_EXTENT(ni1).AND. &
                            & BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2) THEN
                            IF(ni1<ni3) THEN !Special case - put the collapsed node at the end of the line
                              SPECIAL_NODE_COUNT=SPECIAL_NODE_COUNT+1
                              NODES_IN_LINE(MAXIMUM_NODE_EXTENT(ni1))=nn1
                            ELSE
                              NODE_COUNT=NODE_COUNT+1
                              NODES_IN_LINE(NODE_COUNT)=nn1
                            ENDIF
                          ENDIF
                        ELSE
                          !Not collapsed at a xi end
                          IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                            NODE_COUNT=NODE_COUNT+1
                            NODES_IN_LINE(NODE_COUNT)=nn1
                          ENDIF
                        ENDIF
                      ELSE
                        !The current xi must be collapsed
                        IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                          NODE_COUNT=NODE_COUNT+1
                          NODES_IN_LINE(NODE_COUNT)=nn1
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    !The current xi direction, ni1, is not involved in any collapsed (degenerate) planes
                    IF(BASIS%NODE_POSITION_INDEX(nn1,ni2)==nn2.AND.BASIS%NODE_POSITION_INDEX(nn1,ni3)==nn3) THEN
                      NODE_COUNT=NODE_COUNT+1
                      NODES_IN_LINE(NODE_COUNT)=nn1
                    ENDIF
                  ENDIF
                ENDDO !nn1
                IF((NODE_COUNT+SPECIAL_NODE_COUNT)>1) THEN !More than one node so it is a proper line 
                  BASIS%NUMBER_OF_LOCAL_LINES=BASIS%NUMBER_OF_LOCAL_LINES+1
                  BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES)=NODE_COUNT+SPECIAL_NODE_COUNT
                  BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES), &
                    & BASIS%NUMBER_OF_LOCAL_LINES)=NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE( &
                    & BASIS%NUMBER_OF_LOCAL_LINES))
                  DO nnl=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES)
                    DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,BASIS%NUMBER_OF_LOCAL_LINES))
                      IF(BASIS%DERIVATIVE_ORDER_INDEX(nk,BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,BASIS%NUMBER_OF_LOCAL_LINES),ni1)== &
                        & FIRST_PART_DERIV) THEN
                        BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nnl,BASIS%NUMBER_OF_LOCAL_LINES)=nk
                        EXIT
                      ENDIF
                    ENDDO !nk
                  ENDDO !nnl
                  BASIS%LOCAL_LINE_XI_DIRECTION(BASIS%NUMBER_OF_LOCAL_LINES)=ni1
                ENDIF
              ENDDO !nn2
            ENDDO !nn3
         ENDDO !ni1


         !Find the faces
         nn4=1
         ef=0           !element face counter
         DO ni1=1,3
            !Go through all three directions
            ni2=OTHER_XI_DIRECTIONS3(ni1,2,1)
            ni3=OTHER_XI_DIRECTIONS3(ni1,3,1)
            !Pointers for argument list
            np1=>nn1
            np2=>nn2
            np3=>nn3
            np4=>nn4
            IF (ni1==1) THEN
               ARGLIST=INTG_POINTER(np1,np2,np3,np4)
            ELSE IF (ni1==2) THEN
               ARGLIST=INTG_POINTER(np2,np1,np3,np4)
            ELSE
               ARGLIST=INTG_POINTER(np2,np3,np1,np4)
            ENDIF

            nn1=MAXIMUM_NODE_EXTENT(ni1)
            !start with eventually upper face of ni1
            LOCAL_NODE_COUNT=0
            
            IF(BASIS%COLLAPSED_XI(ni1)/=BASIS_COLLAPSED_AT_XI1) THEN
               ef=ef+1  
               DO nn3=1,MAXIMUM_NODE_EXTENT(ni2)
                  DO nn2=1,MAXIMUM_NODE_EXTENT(ni3)
                     LOCAL_NODE_COUNT=LOCAL_NODE_COUNT+1
                     BASIS%NODE_NUMBERS_IN_LOCAL_FACE(LOCAL_NODE_COUNT,ef)= &
                          & BASIS%NODE_POSITION_INDEX_INV(ARGLIST%a1,ARGLIST%a2,ARGLIST%a3,ARGLIST%a4)
                     IF (BASIS%NODE_NUMBERS_IN_LOCAL_FACE(LOCAL_NODE_COUNT,ef)==0) THEN
                        IF (BASIS%COLLAPSED_XI(ni1)==BASIS_XI_COLLAPSED) THEN
                           !set Arglist(ni1-direction)=1
                           nn1=1
                           BASIS%NODE_NUMBERS_IN_LOCAL_FACE(LOCAL_NODE_COUNT,ef)= &
                                & BASIS%NODE_POSITION_INDEX_INV(ARGLIST%a1,ARGLIST%a2,ARGLIST%a3,ARGLIST%a4)
                           nn1=MAXIMUM_NODE_EXTENT(ni1)
                        ELSE
                           LOCAL_NODE_COUNT=LOCAL_NODE_COUNT-1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
               BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(ef)=LOCAL_NODE_COUNT
               BASIS%LOCAL_FACE_XI_DIRECTION(ef)=ni1  
            ENDIF
            
            
            nn1=1
            LOCAL_NODE_COUNT=0
            IF(BASIS%COLLAPSED_XI(ni1)/=BASIS_COLLAPSED_AT_XI0) THEN
               ef=ef+1  
               DO nn3=1,MAXIMUM_NODE_EXTENT(ni2)
                  DO nn2=1,MAXIMUM_NODE_EXTENT(ni3)
                     LOCAL_NODE_COUNT=LOCAL_NODE_COUNT+1
                     BASIS%NODE_NUMBERS_IN_LOCAL_FACE(LOCAL_NODE_COUNT,ef)= &
                          & BASIS%NODE_POSITION_INDEX_INV(ARGLIST%a1,ARGLIST%a2,ARGLIST%a3,ARGLIST%a4)  
                     IF (BASIS%NODE_NUMBERS_IN_LOCAL_FACE(LOCAL_NODE_COUNT,ef)==0) THEN
                        IF (BASIS%COLLAPSED_XI(ni1)==BASIS_XI_COLLAPSED) THEN
                           !set Arglist(ni1-direction)=Max extent
                           nn1=MAXIMUM_NODE_EXTENT(ni1)
                           BASIS%NODE_NUMBERS_IN_LOCAL_FACE(LOCAL_NODE_COUNT,ef)= &
                                & BASIS%NODE_POSITION_INDEX_INV(ARGLIST%a1,ARGLIST%a2,ARGLIST%a3,ARGLIST%a4)
                           nn1=1
                        ELSE
                           LOCAL_NODE_COUNT=LOCAL_NODE_COUNT-1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
               BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(ef)=LOCAL_NODE_COUNT
               BASIS%LOCAL_FACE_XI_DIRECTION(ef)=-ni1  
            ENDIF
         ENDDO
         
         NUMBER_OF_ELEMENT_FACES=ef
         DO ni1=1,3
            DO ef=1,NUMBER_OF_ELEMENT_FACES
               DO nnf=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(ef)
                  DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(nnf,ef))
                     IF(BASIS%DERIVATIVE_ORDER_INDEX(nk,BASIS%NODE_NUMBERS_IN_LOCAL_FACE(nnf,ef),ni1)== &
                          & FIRST_PART_DERIV) THEN
                        BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(nnf,ef)=nk
                     ENDIF
                  ENDDO !nk
               ENDDO !nnf
            ENDDO !ef
         ENDDO !ni1

      CASE DEFAULT
          CALL FLAG_ERROR("Invalid number of xi directions.",ERR,ERROR,*999)
        END SELECT
        
        CALL BASIS_QUADRATURE_CREATE(BASIS,ERR,ERROR,*999)
          
      ELSE
        CALL FLAG_ERROR("Basis is not a Lagrange Hermite tensor product basis.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF

    
    CALL EXITS("BASIS_LHTP_BASIS_CREATE")
    RETURN
999 CALL ERRORS("BASIS_LHTP_BASIS_CREATE",ERR,ERROR)
    CALL EXITS("BASIS_LHTP_BASIS_CREATE")
    RETURN 1
  END SUBROUTINE BASIS_LHTP_BASIS_CREATE

  !
  !
  !================================================================================================================================
  !
  
  !>Evaluates the double precision Lagrange/Hermite/Fourier tensor product basis function for the given BASIS.
  FUNCTION BASIS_LHTP_BASIS_EVALUATE_DP(BASIS,NODE_NUMBER,DERIVATIVE_NUMBER,PARTIAL_DERIV_INDEX,XI,ERR,ERROR)
      
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to evaluate 
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The local node number of the tensor product basis to evaluate
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The local derivative number of the tensor product basis to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX  !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: XI(:) !<The Xi position to evaluate the basis function at
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_LHTP_BASIS_EVALUATE_DP !<On return the evaluated basis funtion.
    !Local variables
    INTEGER(INTG) :: ni,nn
    REAL(DP) :: SUM
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_LHTP_BASIS_EVALUATE_DP",ERR,ERROR,*999)
    
    BASIS_LHTP_BASIS_EVALUATE_DP=1.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      DO ni=1,BASIS%NUMBER_OF_XI
        IF(BASIS%NODE_AT_COLLAPSE(NODE_NUMBER).AND.BASIS%COLLAPSED_XI(ni)==BASIS_XI_COLLAPSED) THEN
          !We are at a collapsed node in the collapsed xi direction. Sum the basis functions in the collapsed xi direction.
          SUM=0.0_DP
          SELECT CASE(BASIS%INTERPOLATION_TYPE(ni))
          CASE(BASIS_LAGRANGE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
              DO nn=1,2
                SUM=SUM+LAGRANGE_LINEAR_EVALUATE(nn,PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
              ENDDO !nn
            CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
              DO nn=1,3
                SUM=SUM+LAGRANGE_QUADRATIC_EVALUATE(nn,PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
              ENDDO !nn
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              DO nn=1,4
                SUM=SUM+LAGRANGE_CUBIC_EVALUATE(nn,PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
              ENDDO !nn
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(ni),"*",ERR,ERROR))// &
                & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(BASIS_HERMITE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              DO nn=1,2
                SUM=SUM+HERMITE_CUBIC_EVALUATE(nn,BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                  & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
              ENDDO !nn
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(ni),"*",ERR,ERROR))// &
                & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            IF(ERR/=0) GOTO 999
          CASE DEFAULT
            LOCAL_ERROR="Interpolation type value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_TYPE(ni),"*",ERR,ERROR))// &
              & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP*SUM
        ELSE
          SELECT CASE(BASIS%INTERPOLATION_TYPE(ni))
          CASE(BASIS_LAGRANGE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & LAGRANGE_LINEAR_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
            CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & LAGRANGE_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & LAGRANGE_CUBIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(ni),"*",ERR,ERROR))// &
                & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(BASIS_HERMITE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_QUADRATIC1_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & HERMITE_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),1,XI(ni),ERR,ERROR)
            CASE(BASIS_QUADRATIC2_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & HERMITE_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),2,XI(ni),ERR,ERROR)
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & HERMITE_CUBIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),ERR,ERROR)
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(ni),"*",ERR,ERROR))// &
                & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            IF(ERR/=0) GOTO 999
          CASE DEFAULT
            LOCAL_ERROR="Interpolation type value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_TYPE(ni),"*",ERR,ERROR))// &
              & " for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ENDDO !ni
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BASIS_LHTP_BASIS_EVALUATE_DP")
    RETURN
999 CALL ERRORS("BASIS_LHTP_BASIS_EVALUATE_DP",ERR,ERROR)
    CALL EXITS("BASIS_LHTP_BASIS_EVALUATE_DP")
    RETURN 
  END FUNCTION BASIS_LHTP_BASIS_EVALUATE_DP

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Lagrange-Hermite tensor product basis family that has already been allocated by BASIS_CREATE_START.
  !>\see BASIS_ROUTINES::BASIS_LHTP_CREATE
  SUBROUTINE BASIS_LHTP_FAMILY_CREATE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,ni,ni2,FACE_XI(2),FACE_XI2(2)
    LOGICAL :: LINE_BASIS_DONE,FACE_BASIS_DONE    
    TYPE(BASIS_TYPE), POINTER :: NEW_SUB_BASIS
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_SUB_BASIS)

    CALL ENTERS("BASIS_LHTP_FAMILY_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      !Create the main (parent) basis
      CALL BASIS_LHTP_BASIS_CREATE(BASIS,ERR,ERROR,*999)
      IF(BASIS%NUMBER_OF_XI>1) THEN
        !Create the line bases as sub-basis types
        ALLOCATE(BASIS%LINE_BASES(BASIS%NUMBER_OF_XI),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis line bases",ERR,ERROR,*999)
        DO ni=1,BASIS%NUMBER_OF_XI
          LINE_BASIS_DONE=.FALSE.
          NULLIFY(NEW_SUB_BASIS)
          DO ni2=1,ni-1
            IF(BASIS%INTERPOLATION_XI(ni2)==BASIS%INTERPOLATION_XI(ni).AND. &
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni2)==BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)) THEN
              LINE_BASIS_DONE=.TRUE.
              EXIT
            ENDIF
          ENDDO !ni2
          IF(LINE_BASIS_DONE) THEN
            BASIS%LINE_BASES(ni)%PTR=>BASIS%LINE_BASES(ni2)%PTR
          ELSE
            !Create the new sub-basis
            CALL BASIS_SUB_BASIS_CREATE(BASIS,1,(/ni/),NEW_SUB_BASIS,ERR,ERROR,*999)
            !Fill in the basis information
            CALL BASIS_LHTP_BASIS_CREATE(NEW_SUB_BASIS,ERR,ERROR,*999)
            BASIS%LINE_BASES(ni)%PTR=>NEW_SUB_BASIS
          ENDIF
        ENDDO !ni
        IF(BASIS%NUMBER_OF_XI>2) THEN
          !Set up face basis functions
          ALLOCATE(BASIS%FACE_BASES(BASIS%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis face bases",ERR,ERROR,*999)
          DO ni=1,BASIS%NUMBER_OF_XI
            !Determine the face xi directions that lie in this xi direction
            FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
            FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
            FACE_BASIS_DONE=.FALSE.
            NULLIFY(NEW_SUB_BASIS)
            DO ni2=1,ni-1
              FACE_XI2(1)=OTHER_XI_DIRECTIONS3(ni2,2,1)
              FACE_XI2(2)=OTHER_XI_DIRECTIONS3(ni2,3,1)
              !Going to disable the test below, as it results in error in collapsed elements and doesn't save much time
!               IF(BASIS%INTERPOLATION_XI(FACE_XI2(1))==BASIS%INTERPOLATION_XI(FACE_XI(1)).AND. &
!                 & BASIS%INTERPOLATION_XI(FACE_XI2(2))==BASIS%INTERPOLATION_XI(FACE_XI(2)).AND. &
!                 & BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI2(1))==BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(1)).AND. &
!                 & BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI2(2))==BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(1))) THEN
!                 FACE_BASIS_DONE=.TRUE.
!                 EXIT
!               ENDIF
            ENDDO !ni2
            IF(FACE_BASIS_DONE) THEN
              BASIS%FACE_BASES(ni)%PTR=>BASIS%FACE_BASES(ni2)%PTR
            ELSE
              !Create the new sub-basis
              CALL BASIS_SUB_BASIS_CREATE(BASIS,2,(/FACE_XI(1),FACE_XI(2)/),NEW_SUB_BASIS,ERR,ERROR,*999)
              !Fill in the basis information
              CALL BASIS_LHTP_BASIS_CREATE(NEW_SUB_BASIS,ERR,ERROR,*999)
              NEW_SUB_BASIS%LINE_BASES(1)%PTR=>BASIS%LINE_BASES(FACE_XI(1))%PTR
              NEW_SUB_BASIS%LINE_BASES(2)%PTR=>BASIS%LINE_BASES(FACE_XI(2))%PTR
              BASIS%FACE_BASES(ni)%PTR=>NEW_SUB_BASIS
            ENDIF            
          ENDDO !ni
        ELSE
          ALLOCATE(BASIS%FACE_BASES(1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis face bases",ERR,ERROR,*999)
          BASIS%FACE_BASES(1)%PTR=>BASIS
        ENDIF
      ELSE
        ALLOCATE(BASIS%LINE_BASES(1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis line bases",ERR,ERROR,*999)
        BASIS%LINE_BASES(1)%PTR=>BASIS
        NULLIFY(BASIS%FACE_BASES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_LHTP_FAMILY_CREATE")
    RETURN
999 IF(ASSOCIATED(NEW_SUB_BASIS)) CALL BASIS_FAMILY_DESTROY(NEW_SUB_BASIS%USER_NUMBER,NEW_SUB_BASIS%FAMILY_NUMBER, &
      & DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("BASIS_LHTP_FAMILY_CREATE",ERR,ERROR)
    CALL EXITS("BASIS_LHTP_FAMILY_CREATE")
    RETURN 1
  END SUBROUTINE BASIS_LHTP_FAMILY_CREATE

  !
  !================================================================================================================================
  !

  !>Calculates the xi location of a local node in a basis.
  SUBROUTINE BASIS_LOCAL_NODE_XI_CALCULATE(BASIS,LOCAL_NODE_NUMBER,XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to calculate the xi for
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to calculate the xi for
    REAL(DP), INTENT(OUT) :: XI(:) !<XI(ni). On return, the xi position of the local node in the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: xi_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_LOCAL_NODE_XI_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=BASIS%NUMBER_OF_NODES) THEN
          IF(SIZE(XI,1)>=BASIS%NUMBER_OF_XI) THEN
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO xi_idx=1,BASIS%NUMBER_OF_XI
                XI(xi_idx)=REAL(BASIS%NODE_POSITION_INDEX(LOCAL_NODE_NUMBER,xi_idx)-1,DP)/ &
                  & REAL(BASIS%NUMBER_OF_NODES_XIC(xi_idx)-1,DP)
              ENDDO !xi_idx
            CASE(BASIS_SIMPLEX_TYPE)
              DO xi_idx=1,BASIS%NUMBER_OF_XI
                XI(xi_idx)=REAL(BASIS%NUMBER_OF_NODES_XIC(xi_idx)-BASIS%NODE_POSITION_INDEX(LOCAL_NODE_NUMBER,xi_idx),DP)/ &
                  & REAL(BASIS%NUMBER_OF_NODES_XIC(xi_idx)-1,DP)
              ENDDO !xi_idx
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The size of the specified xic array of "//TRIM(NUMBER_TO_VSTRING(SIZE(XI,1),"*",ERR,ERROR))// &
              & " is invalid. The size of the xi array must be >= "// &
              & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."            
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified local node number of "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The local node number must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Basis has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_LOCAL_NODE_XI_CALCULATE")
    RETURN
999 CALL ERRORS("BASIS_LOCAL_NODE_XI_CALCULATE",ERR,ERROR)
    CALL EXITS("BASIS_LOCAL_NODE_XI_CALCULATE")
    RETURN 1
  END SUBROUTINE BASIS_LOCAL_NODE_XI_CALCULATE

  !
  !================================================================================================================================
  !

  !>Returns the number of local nodes in the specified basis \see OPENCMISS::CMISSBasisNumberOfLocalNodesGet
  SUBROUTINE BASIS_NUMBER_OF_LOCAL_NODES_GET(BASIS,NUMBER_OF_LOCAL_NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to get the number of nodes
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_LOCAL_NODES !<On return, the number of local nodes in the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BASIS_NUMBER_OF_LOCAL_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      NUMBER_OF_LOCAL_NODES=BASIS%NUMBER_OF_NODES
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_NUMBER_OF_LOCAL_NODES_GET")
    RETURN
999 CALL ERRORS("BASIS_NUMBER_OF_LOCAL_NODES_GET",ERR,ERROR)
    CALL EXITS("BASIS_NUMBER_OF_LOCAL_NODES_GET")
    RETURN 1
  END SUBROUTINE BASIS_NUMBER_OF_LOCAL_NODES_GET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of xi directions for a basis. \see OPENCMISS::CMISSBasisNumberOfXiGet
  SUBROUTINE BASIS_NUMBER_OF_XI_GET(BASIS,NUMBER_OF_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function to change
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_XI !<On return, the number of Xi directions to get.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BASIS_NUMBER_OF_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        NUMBER_OF_XI=BASIS%NUMBER_OF_XI
      ELSE
        CALL FLAG_ERROR("Basis has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_NUMBER_OF_XI_GET")
    RETURN
999 CALL ERRORS("BASIS_NUMBER_OF_XI_GET",ERR,ERROR)
    CALL EXITS("BASIS_NUMBER_OF_XI_GET")
    RETURN
  END SUBROUTINE BASIS_NUMBER_OF_XI_GET 

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of xi directions where the basis is identified by user number.
  SUBROUTINE BASIS_NUMBER_OF_XI_SET_NUMBER(USER_NUMBER,NUMBER_OF_XI,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI !<The number of Xi directions to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_NUMBER_OF_XI_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_NUMBER_OF_XI_SET(BASIS,NUMBER_OF_XI,ERR,ERROR,*999)

    CALL EXITS("BASIS_NUMBER_OF_XI_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_NUMBER_OF_XI_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_NUMBER_OF_XI_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_NUMBER_OF_XI_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of xi directions for a basis identified by a pointer. \see OPENCMISS::CMISSBasisNumberOfXiSet
  SUBROUTINE BASIS_NUMBER_OF_XI_SET_PTR(BASIS,NUMBER_OF_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function to change
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI !<The number of Xi directions to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: OLD_INTERPOLATION_XI(3),OLD_NUMBER_OF_GAUSS_XI(3),OLD_COLLAPSED_XI(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_NUMBER_OF_XI_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(BASIS%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          IF(NUMBER_OF_XI>0.AND.NUMBER_OF_XI<4) THEN
            IF(BASIS%NUMBER_OF_XI/=NUMBER_OF_XI) THEN
              !Reallocate the basis information arrays that depend on the number of xi directions
              OLD_INTERPOLATION_XI=BASIS%INTERPOLATION_XI
              OLD_COLLAPSED_XI=BASIS%COLLAPSED_XI
              DEALLOCATE(BASIS%INTERPOLATION_XI)
              DEALLOCATE(BASIS%COLLAPSED_XI)
              ALLOCATE(BASIS%INTERPOLATION_XI(NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type",ERR,ERROR,*999)
              ALLOCATE(BASIS%COLLAPSED_XI(NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate collapsed xi",ERR,ERROR,*999)
              IF(NUMBER_OF_XI>BASIS%NUMBER_OF_XI) THEN
                BASIS%INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)=OLD_INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)
                BASIS%INTERPOLATION_XI(BASIS%NUMBER_OF_XI+1:NUMBER_OF_XI)=OLD_INTERPOLATION_XI(1)
                BASIS%COLLAPSED_XI(1:BASIS%NUMBER_OF_XI)=OLD_COLLAPSED_XI(1:BASIS%NUMBER_OF_XI)
                BASIS%COLLAPSED_XI(BASIS%NUMBER_OF_XI+1:NUMBER_OF_XI)=OLD_COLLAPSED_XI(1)
              ELSE
                BASIS%INTERPOLATION_XI(1:NUMBER_OF_XI)=OLD_INTERPOLATION_XI(1:NUMBER_OF_XI)
                BASIS%COLLAPSED_XI(1:NUMBER_OF_XI)=OLD_COLLAPSED_XI(1:NUMBER_OF_XI)
              ENDIF
              
              IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
                OLD_NUMBER_OF_GAUSS_XI=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI
                DEALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)
                ALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(NUMBER_OF_XI),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of Gauss xi",ERR,ERROR,*999)
                IF(NUMBER_OF_XI>BASIS%NUMBER_OF_XI) THEN
                  BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI)=OLD_NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI)
                  BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(BASIS%NUMBER_OF_XI+1:NUMBER_OF_XI)=OLD_NUMBER_OF_GAUSS_XI(1)
                ELSE
                  BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:NUMBER_OF_XI)=OLD_NUMBER_OF_GAUSS_XI(1:NUMBER_OF_XI)
                ENDIF
              ENDIF
              BASIS%NUMBER_OF_XI=NUMBER_OF_XI
            ENDIF
          ELSE
            LOCAL_ERROR="Invalid number of xi directions specified ("// &
              & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_XI,"*",ERR,ERROR))// &
              & ") for a Lagrange-Hermite basis. You must specify between 1 and 3 xi directions"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(BASIS_SIMPLEX_TYPE)
          IF(NUMBER_OF_XI>1.AND.NUMBER_OF_XI<4) THEN
            IF(BASIS%NUMBER_OF_XI/=NUMBER_OF_XI) THEN
              !Reallocate the basis information arrays that depend on the number of xi directions
              OLD_INTERPOLATION_XI=BASIS%INTERPOLATION_XI
              OLD_COLLAPSED_XI=BASIS%COLLAPSED_XI
              DEALLOCATE(BASIS%INTERPOLATION_XI)
              DEALLOCATE(BASIS%COLLAPSED_XI)
              ALLOCATE(BASIS%INTERPOLATION_XI(NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type",ERR,ERROR,*999)
              ALLOCATE(BASIS%COLLAPSED_XI(NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate collapsed xi.",ERR,ERROR,*999)
              IF(NUMBER_OF_XI>BASIS%NUMBER_OF_XI) THEN
                BASIS%INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)=OLD_INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)
                BASIS%INTERPOLATION_XI(BASIS%NUMBER_OF_XI+1:NUMBER_OF_XI)=OLD_INTERPOLATION_XI(1)
                BASIS%COLLAPSED_XI(1:BASIS%NUMBER_OF_XI)=OLD_COLLAPSED_XI(1:BASIS%NUMBER_OF_XI)
                BASIS%COLLAPSED_XI(BASIS%NUMBER_OF_XI+1:NUMBER_OF_XI)=OLD_COLLAPSED_XI(1)
              ELSE
                BASIS%INTERPOLATION_XI(1:NUMBER_OF_XI)=OLD_INTERPOLATION_XI(1:NUMBER_OF_XI)
                BASIS%COLLAPSED_XI(1:NUMBER_OF_XI)=OLD_COLLAPSED_XI(1:NUMBER_OF_XI)
              ENDIF              
              IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
                OLD_NUMBER_OF_GAUSS_XI=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI
                DEALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)
                ALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(NUMBER_OF_XI),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of Gauss xi",ERR,ERROR,*999)
                IF(NUMBER_OF_XI>BASIS%NUMBER_OF_XI) THEN
                  BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI)=OLD_NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI)
                  BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(BASIS%NUMBER_OF_XI+1:NUMBER_OF_XI)=OLD_NUMBER_OF_GAUSS_XI(1)
                ELSE
                  BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:NUMBER_OF_XI)=OLD_NUMBER_OF_GAUSS_XI(1:NUMBER_OF_XI)
                ENDIF
              ENDIF
              BASIS%NUMBER_OF_XI=NUMBER_OF_XI
            ENDIF
          ELSE
            LOCAL_ERROR="Invalid number of xi directions specified ("// &
              & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_XI,"*",ERR,ERROR))// &
              & ") for a simplex basis. You must specify between 2 and 3 xi directions"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Basis type invalid or not implemented",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_NUMBER_OF_XI_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_NUMBER_OF_XI_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_NUMBER_OF_XI_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_NUMBER_OF_XI_SET_PTR 
  
  !
  !================================================================================================================================
  !

  !>Creates the quadrature and quadrature schemes on a basis.
  SUBROUTINE BASIS_QUADRATURE_CREATE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: scheme_idx,i,j,k,MAX_NUM_GAUSS,ng,ni,nk,nn,ns,nu,NUM_GAUSS_1,NUM_GAUSS_2,NUM_GAUSS_3
    REAL(DP) :: XI(3),GSX(4,20),GSW(20)
    REAL(DP), ALLOCATABLE :: POSITIONS(:,:),POSITIONS_MATRIX(:,:,:,:),WEIGHTS(:,:)
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: NEW_SCHEME,SCHEME
    TYPE(QUADRATURE_SCHEME_PTR_TYPE), POINTER :: NEW_SCHEMES(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: MAX_NUM_FACE_GAUSS,face_idx,NORMAL,FACE_XI(2)

    NULLIFY(NEW_SCHEME)
    NULLIFY(NEW_SCHEMES)

    CALL ENTERS("BASIS_QUADRATURE_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ASSOCIATED(BASIS%QUADRATURE%SCHEMES)) THEN
        LOCAL_ERROR="The quadrature schemes on basis number "//TRIM(NUMBER_TO_VSTRING(BASIS%USER_NUMBER,"*",ERR,ERROR))// &
          & " are already associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ELSE
!!TODO: \todo Sort this properly by having a create values cache.
        !Reset the basis quadrature - 
        !CALL BASIS_QUADRATURE_FINALISE(BASIS,ERR,ERROR,*999)      ! Kumar - I think this is not correct as it 
        !Initialise the basis quadrature                           !         resets the quadrature scheme already set.
        !CALL BASIS_QUADRATURE_INITIALISE(BASIS,ERR,ERROR,*999)    ! 
        SELECT CASE(BASIS%QUADRATURE%TYPE)
        CASE(BASIS_GAUSS_LEGENDRE_QUADRATURE)            
          !Allocate one scheme and add it to the list of schemes
          ALLOCATE(NEW_SCHEME,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new quadrature scheme",ERR,ERROR,*999)
          NEW_SCHEME%QUADRATURE=>BASIS%QUADRATURE
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=1
          ALLOCATE(NEW_SCHEMES(BASIS%QUADRATURE%NUMBER_OF_SCHEMES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new quadratures scheme",ERR,ERROR,*999)
          NEW_SCHEMES(1)%PTR=>NEW_SCHEME
          BASIS%QUADRATURE%SCHEMES=>NEW_SCHEMES
          !Set up the quadrature scheme map
          ALLOCATE(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate quadrature scheme map",ERR,ERROR,*999)
          DO scheme_idx=1,BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES
            NULLIFY(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(scheme_idx)%PTR)
          ENDDO !scheme_idx
          BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR=>NEW_SCHEME
          !Set up the gauss point arrays
          NEW_SCHEME%NUMBER_OF_GAUSS=1
          MAX_NUM_GAUSS=-1
          DO ni=1,BASIS%NUMBER_OF_XI
            NEW_SCHEME%NUMBER_OF_GAUSS=NEW_SCHEME%NUMBER_OF_GAUSS*BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)
            IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)>MAX_NUM_GAUSS) MAX_NUM_GAUSS=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)
          ENDDO !ni
          ALLOCATE(NEW_SCHEME%GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Gauss positions",ERR,ERROR,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_WEIGHTS(NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Gauss weights",ERR,ERROR,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
            & NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Gauss basis functions",ERR,ERROR,*999)
          ALLOCATE(WEIGHTS(MAX_NUM_GAUSS,3),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate weights",ERR,ERROR,*999)
          ALLOCATE(POSITIONS(MAX_NUM_GAUSS,3),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate positions",ERR,ERROR,*999)
          ALLOCATE(POSITIONS_MATRIX(MAX_NUM_GAUSS,MAX_NUM_GAUSS,MAX_NUM_GAUSS,3),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate positions matrix",ERR,ERROR,*999)
          WEIGHTS=1.0_DP
          POSITIONS=0.0_DP
          POSITIONS_MATRIX=0.0_DP
          DO ni=1,BASIS%NUMBER_OF_XI
            CALL GAUSS_LEGENDRE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),0.0_DP,1.0_DP, &
              & POSITIONS(1:BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),ni), &
              & WEIGHTS(1:BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),ni),ERR,ERROR,*999)
          ENDDO !ni
          SELECT CASE(BASIS%NUMBER_OF_XI)
          CASE(1)
            NUM_GAUSS_1=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1)
            NUM_GAUSS_2=1
            NUM_GAUSS_3=1
          CASE(2)
            NUM_GAUSS_1=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1)
            NUM_GAUSS_2=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(2)
            NUM_GAUSS_3=1
          CASE(3)
            NUM_GAUSS_1=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1)
            NUM_GAUSS_2=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(2)
            NUM_GAUSS_3=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(3)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid number of xi directions",ERR,ERROR,*999)
          END SELECT
          DO k=1,NUM_GAUSS_3
            DO j=1,NUM_GAUSS_2
              DO i=1,NUM_GAUSS_1
                POSITIONS_MATRIX(i,j,k,1)=POSITIONS(i,1)
                POSITIONS_MATRIX(i,j,k,2)=POSITIONS(j,2)
                POSITIONS_MATRIX(i,j,k,3)=POSITIONS(k,3)
                XI(1:BASIS%NUMBER_OF_XI)=POSITIONS_MATRIX(i,j,k,1:BASIS%NUMBER_OF_XI)
                ng=i+(j-1+(k-1)*NUM_GAUSS_2)*NUM_GAUSS_1
                NEW_SCHEME%GAUSS_WEIGHTS(ng)=WEIGHTS(i,1)*WEIGHTS(j,2)*WEIGHTS(k,3)
                NEW_SCHEME%GAUSS_POSITIONS(1:BASIS%NUMBER_OF_XI_COORDINATES,ng)=XI(1:BASIS%NUMBER_OF_XI_COORDINATES)
                ns=0
                DO nn=1,BASIS%NUMBER_OF_NODES
                  DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                    ns=ns+1
                    DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                      SELECT CASE(BASIS%TYPE)
                      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                        NEW_SCHEME%GAUSS_BASIS_FNS(ns,nu,ng)=BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,nu,XI,ERR,ERROR)
                        IF(ERR/=0) GOTO 999                        
                      CASE DEFAULT
                        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                      END SELECT
                    ENDDO !nu
                  ENDDO !nk
                ENDDO !nn
              ENDDO !i
            ENDDO !j
          ENDDO !k
          !Create face quadrature scheme, if requested
          IF(BASIS%QUADRATURE%EVALUATE_FACE_GAUSS) THEN
            IF(BASIS%NUMBER_OF_XI==3) THEN
              !Find maximum number of face gauss points and allocate the arrays
              MAX_NUM_FACE_GAUSS=PRODUCT(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI))
              MAX_NUM_FACE_GAUSS=MAX_NUM_FACE_GAUSS/MINVAL(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI))
              ALLOCATE(NEW_SCHEME%NUMBER_OF_FACE_GAUSS(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of face gauss",ERR,ERROR,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,MAX_NUM_FACE_GAUSS, &
                & BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face Gauss positions",ERR,ERROR,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_WEIGHTS(MAX_NUM_FACE_GAUSS,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face Gauss weights",ERR,ERROR,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
                & MAX_NUM_FACE_GAUSS,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face Gauss basis function values array",ERR,ERROR,*999)
              !Zero them out just to be safe
              NEW_SCHEME%FACE_GAUSS_POSITIONS=0.0_DP
              NEW_SCHEME%FACE_GAUSS_WEIGHTS=0.0_DP
              NEW_SCHEME%FACE_GAUSS_BASIS_FNS=0.0_DP
              !Populate face_gauss_positions, weights, basis_fn
              DO face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                !What's the normal?
                NORMAL=BASIS%LOCAL_FACE_XI_DIRECTION(face_idx)
                IF(NORMAL<0_INTG) THEN
                  XI(ABS(NORMAL))=0.0_DP
                ELSE
                  XI(ABS(NORMAL))=1.0_DP
                ENDIF
                NORMAL=ABS(NORMAL)
                FACE_XI=[OTHER_XI_DIRECTIONS3(NORMAL,2,1), OTHER_XI_DIRECTIONS3(NORMAL,3,1)]
                !How many gauss points are in this face?
                NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx)=PRODUCT(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI))
                ng=0_INTG
                DO j=1,BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(2))
                  XI(FACE_XI(2))=POSITIONS(j,FACE_XI(2))
                  DO i=1,BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(1))
                    XI(FACE_XI(1))=POSITIONS(i,FACE_XI(1))
                    ng=ng+1_INTG
                    !Gauss point xi and weights first
                    NEW_SCHEME%FACE_GAUSS_WEIGHTS(ng,face_idx)=WEIGHTS(i,FACE_XI(1))*WEIGHTS(j,FACE_XI(2))
                    NEW_SCHEME%FACE_GAUSS_POSITIONS(1:3,ng,face_idx)=XI(1:3)
                    !Evaluate basis fn values at the Gauss points now
                    ns=0
                    DO nn=1,BASIS%NUMBER_OF_NODES
                      DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                        ns=ns+1
                        DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                          SELECT CASE(BASIS%TYPE)
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                            NEW_SCHEME%FACE_GAUSS_BASIS_FNS(ns,nu,ng,face_idx)= &
                              & BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,nu,XI,ERR,ERROR)
                            IF(ERR/=0) GOTO 999                        
                          CASE DEFAULT
                            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                          END SELECT
                        ENDDO !nu
                      ENDDO !nk
                    ENDDO !nn

                  ENDDO !i
                ENDDO !j
              ENDDO !face_idx
            ELSE
              CALL FLAG_ERROR("Cannot evaluate face quadrature schemes for a non three dimensional element.",ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Clean up
          DEALLOCATE(WEIGHTS)
          DEALLOCATE(POSITIONS)
          DEALLOCATE(POSITIONS_MATRIX)
        CASE(BASIS_GAUSS_LAGUERRE_QUADRATURE)
          CALL FLAG_ERROR("Gauss Laguerre quadrature type not implemented",ERR,ERROR,*999)
        CASE(BASIS_GUASS_HERMITE_QUADRATURE)
          CALL FLAG_ERROR("Gauss Hermite quadrature type not implemented",ERR,ERROR,*999)
        CASE(BASIS_GAUSS_SIMPLEX_QUADRATURE)
          !Allocate one scheme and add it to the list of schemes
          ALLOCATE(NEW_SCHEME,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new quadrature scheme",ERR,ERROR,*999)
          NEW_SCHEME%QUADRATURE=>BASIS%QUADRATURE
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=1
          ALLOCATE(NEW_SCHEMES(BASIS%QUADRATURE%NUMBER_OF_SCHEMES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new quadratures scheme",ERR,ERROR,*999)
          NEW_SCHEMES(1)%PTR=>NEW_SCHEME
          BASIS%QUADRATURE%SCHEMES=>NEW_SCHEMES
          !Set up the quadrature scheme map
          ALLOCATE(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate quadrature scheme map",ERR,ERROR,*999)
          DO scheme_idx=1,BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES
            NULLIFY(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(scheme_idx)%PTR)
          ENDDO !scheme_idx
          BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR=>NEW_SCHEME
          !Set up the gauss point arrays
          CALL GAUSS_SIMPLEX(BASIS%QUADRATURE%GAUSS_ORDER,BASIS%NUMBER_OF_XI_COORDINATES,NEW_SCHEME%NUMBER_OF_GAUSS,GSX,GSW, &
            & ERR,ERROR,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Gauss positions",ERR,ERROR,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_WEIGHTS(NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Gauss weights",ERR,ERROR,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
            & NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Gauss basis functions",ERR,ERROR,*999)
          NEW_SCHEME%GAUSS_POSITIONS=GSX(1:BASIS%NUMBER_OF_XI_COORDINATES,1:NEW_SCHEME%NUMBER_OF_GAUSS)
          NEW_SCHEME%GAUSS_WEIGHTS=GSW(1:NEW_SCHEME%NUMBER_OF_GAUSS)
          DO ng=1,NEW_SCHEME%NUMBER_OF_GAUSS
            ns=0
            DO nn=1,BASIS%NUMBER_OF_NODES
              DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                ns=ns+1
                DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                  SELECT CASE(BASIS%TYPE)
                  CASE(BASIS_SIMPLEX_TYPE)
                    NEW_SCHEME%GAUSS_BASIS_FNS(ns,nu,ng)= &
                      & BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,nu,NEW_SCHEME%GAUSS_POSITIONS(1:BASIS%NUMBER_OF_XI_COORDINATES,ng), &
                      & ERR,ERROR)
                    IF(ERR/=0) GOTO 999                        
                  CASE DEFAULT
                    CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                  END SELECT
                ENDDO !nu
              ENDDO !nk
            ENDDO !nn
          ENDDO !ng
        CASE DEFAULT
          LOCAL_ERROR="Quadrature type "//TRIM(NUMBER_TO_VSTRING(BASIS%QUADRATURE%TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Quadrature type = ",BASIS%QUADRATURE%TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI,3,3,BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI, &
        & '("  Number of gauss points(ni):",3(X,I2))','(22X,3(X,I2))',ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of quadrature schemes = ",BASIS%QUADRATURE%NUMBER_OF_SCHEMES, &
        & ERR,ERROR,*999)
      DO scheme_idx=1,BASIS%QUADRATURE%NUMBER_OF_SCHEMES
        SCHEME=>BASIS%QUADRATURE%SCHEMES(scheme_idx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Scheme = ",scheme_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of gauss points = ",SCHEME%NUMBER_OF_GAUSS, &
          & ERR,ERROR,*999)
        IF(DIAGNOSTICS2) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Gauss point positions and weights:",ERR,ERROR,*999)
          DO ng=1,SCHEME%NUMBER_OF_GAUSS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss point = ",ng,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI_COORDINATES,3,3,SCHEME%GAUSS_POSITIONS(:,ng), &
              & '("          POSITION(ni)   :",3(X,F12.4))','(26X,3(X,F12.4))',ERR,ERROR,*999)
            CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          WEIGHT         : ",SCHEME%GAUSS_WEIGHTS(ng), &
              & "(F12.4)",ERR,ERROR,*999)
          ENDDO !ng          
        ENDIF
        IF(DIAGNOSTICS3) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Basis functions evaluated at Gauss points:",ERR,ERROR,*999)
          DO ng=1,SCHEME%NUMBER_OF_GAUSS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss point = ",ng,ERR,ERROR,*999)
            DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Partial derivative number = ",nu,ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS,4,4, &
                & SCHEME%GAUSS_BASIS_FNS(:,nu,ng),'("          BASIS FNS(ns)  :",4(X,F12.4))','(26X,4(X,F12.4))',ERR,ERROR,*999)
            ENDDO !nu
          ENDDO !ng
        ENDIF
      ENDDO !scheme_idx
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_CREATE")
    RETURN
999 IF(ASSOCIATED(NEW_SCHEME)) THEN
      IF(ALLOCATED(NEW_SCHEME%GAUSS_POSITIONS)) DEALLOCATE(NEW_SCHEME%GAUSS_POSITIONS)
      IF(ALLOCATED(NEW_SCHEME%GAUSS_WEIGHTS)) DEALLOCATE(NEW_SCHEME%GAUSS_WEIGHTS)
      IF(ALLOCATED(NEW_SCHEME%GAUSS_BASIS_FNS)) DEALLOCATE(NEW_SCHEME%GAUSS_BASIS_FNS)
      DEALLOCATE(NEW_SCHEME)     
    ENDIF
    IF(ALLOCATED(WEIGHTS)) DEALLOCATE(WEIGHTS)
    IF(ALLOCATED(POSITIONS)) DEALLOCATE(POSITIONS)
    IF(ALLOCATED(POSITIONS_MATRIX)) DEALLOCATE(POSITIONS_MATRIX)
    IF(ASSOCIATED(NEW_SCHEMES)) DEALLOCATE(NEW_SCHEMES)
    NULLIFY(BASIS%QUADRATURE%SCHEMES)
998 CALL ERRORS("BASIS_QUADRATURE_CREATE",ERR,ERROR)    
    CALL EXITS("BASIS_QUADRATURE_CREATE")
    RETURN 1
   
  END SUBROUTINE BASIS_QUADRATURE_CREATE
        
  !
  !================================================================================================================================
  !
  
  !>Destroys a quadrature on a given basis and deallocates all memory. \todo fix all this basis/quadrature into standard form.
  SUBROUTINE BASIS_QUADRATURE_DESTROY(QUADRATURE,ERR,ERROR,*)

    !Argument variables
    TYPE(QUADRATURE_TYPE), POINTER :: QUADRATURE !<A pointer to the quadrature
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASIS_QUADRATURE_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(QUADRATURE)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)     
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_DESTROY")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_DESTROY",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_DESTROY")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_DESTROY

  !
  !================================================================================================================================
  !
    
  !>Finalises a quadrature on a given basis and deallocates all memory
  SUBROUTINE BASIS_QUADRATURE_FINALISE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: scheme_idx
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: SCHEME

    CALL ENTERS("BASIS_QUADRATURE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
        DO scheme_idx=1,BASIS%QUADRATURE%NUMBER_OF_SCHEMES
          SCHEME=>BASIS%QUADRATURE%SCHEMES(scheme_idx)%PTR
          !Destroy all scheme components
          IF(ALLOCATED(SCHEME%GAUSS_POSITIONS)) DEALLOCATE(SCHEME%GAUSS_POSITIONS)
          IF(ALLOCATED(SCHEME%GAUSS_WEIGHTS)) DEALLOCATE(SCHEME%GAUSS_WEIGHTS)
          IF(ALLOCATED(SCHEME%GAUSS_BASIS_FNS)) DEALLOCATE(SCHEME%GAUSS_BASIS_FNS)
        ENDDO !scheme_idx
        IF(ASSOCIATED(BASIS%QUADRATURE%SCHEMES)) DEALLOCATE(BASIS%QUADRATURE%SCHEMES)
        BASIS%QUADRATURE%NUMBER_OF_SCHEMES=0
        IF(ALLOCATED(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)) DEALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)
        NULLIFY(BASIS%QUADRATURE%BASIS)
        BASIS%QUADRATURE%TYPE=-1
      ELSE
        CALL FLAG_ERROR("Basis quadrature basis is not associated",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_FINALISE")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_FINALISE",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_FINALISE")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a quadrature on the given basis.
  SUBROUTINE BASIS_QUADRATURE_INITIALISE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BASIS_QUADRATURE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
        LOCAL_ERROR="Basis number "//TRIM(NUMBER_TO_VSTRING(BASIS%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a quadrature associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ELSE
        SELECT CASE(BASIS%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          !Set up a default Gauss Legendre quadrature 
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=0
          NULLIFY(BASIS%QUADRATURE%SCHEMES)
          BASIS%QUADRATURE%BASIS=>BASIS        
          BASIS%QUADRATURE%TYPE=BASIS_GAUSS_LEGENDRE_QUADRATURE
          !Set up a default number of Gauss points appropriate for the given interpolation order in each direction.
          ALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(BASIS%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of Gauss in each xi direction",ERR,ERROR,*999)
          DO ni=1,BASIS%NUMBER_OF_XI
            SELECT CASE(BASIS%INTERPOLATION_XI(ni))
            CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=2
            CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=3
            CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=4
            CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=3
            CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=4
            CASE DEFAULT
              LOCAL_ERROR="Interpolation xi value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_XI(ni),"*",ERR,ERROR))// &
                & " in xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDDO !ni
        CASE(BASIS_SIMPLEX_TYPE)
         !Set up a default quadrature 
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=0
          NULLIFY(BASIS%QUADRATURE%SCHEMES)
          BASIS%QUADRATURE%BASIS=>BASIS        
          BASIS%QUADRATURE%TYPE=BASIS_GAUSS_SIMPLEX_QUADRATURE
          !Set up a default order appropriate for the given interpolation.
          ALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(BASIS%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of Gauss in each xi direction",ERR,ERROR,*999)
!!TODO: \todo Set these to something more meaningfull!
          SELECT CASE(BASIS%INTERPOLATION_XI(1))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%QUADRATURE%GAUSS_ORDER=2
            CASE(2)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(3)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE DEFAULT
              LOCAL_ERROR="The number of xi directions ("//TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))// &
                & ") is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(2)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(3)
              BASIS%QUADRATURE%GAUSS_ORDER=5
            CASE DEFAULT
              LOCAL_ERROR="The number of xi directions ("//TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))// &
                & ") is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(2)
              BASIS%QUADRATURE%GAUSS_ORDER=5
            CASE(3)
              BASIS%QUADRATURE%GAUSS_ORDER=5
            CASE DEFAULT
              LOCAL_ERROR="The number of xi directions ("//TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))// &
                & ") is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Interpolation xi value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_XI(1),"*",ERR,ERROR))// &
              & " in xi direction 1 is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI=BASIS%QUADRATURE%GAUSS_ORDER
        CASE DEFAULT
          LOCAL_ERROR="Basis type value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_XI(ni),"*",ERR,ERROR))// &
            & " is invalid or not implemented"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_INITIALISE")
    RETURN
999 IF(ALLOCATED(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)) DEALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)
998 CALL ERRORS("BASIS_QUADRATURE_INITIALISE",ERR,ERROR)    
    CALL EXITS("BASIS_QUADRATURE_INITIALISE")
    RETURN 1
    
  END SUBROUTINE BASIS_QUADRATURE_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Get the number of Gauss points in each xi direction on a basis quadrature identified by a pointer. \see OPENCMISS::CMISSBasisQuadratureNumberOfGaussXiGet
  SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET(BASIS,QUADRATURE_NUMBER_OF_GAUSS_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: QUADRATURE_NUMBER_OF_GAUSS_XI(:) !<On return, the number of Gauss in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          IF(SIZE(QUADRATURE_NUMBER_OF_GAUSS_XI,1)>=SIZE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI,1)) THEN
            QUADRATURE_NUMBER_OF_GAUSS_XI=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI
          ELSE
            LOCAL_ERROR="The size of QUADRATURE_NUMBER_OF_GAUSS_XI is too small. The supplied size is "// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(QUADRATURE_NUMBER_OF_GAUSS_XI,1),"*",ERR,ERROR))//" and it needs to be >= "// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI,1),"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Quadrature basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Basis has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET")
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET
    
  !
  !================================================================================================================================
  !

  !>Sets/changes the number of gauss in each xi direction for a basis quadrature identified by a user number.
  SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER(USER_NUMBER,NUMBER_GAUSS_XI,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis
    INTEGER(INTG), INTENT(IN) :: NUMBER_GAUSS_XI(:) !<The number of Gauss points in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(BASIS,NUMBER_GAUSS_XI,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of Gauss points in each xi direction on a basis quadrature identified by a pointer. \see OPENCMISS::CMISSBasisQuadratureNumberOfGaussSet
  SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR(BASIS,NUMBER_OF_GAUSS_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_GAUSS_XI(:) !<The number of Gauss in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_WARNING
    
    CALL ENTERS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN          
          IF(SIZE(NUMBER_OF_GAUSS_XI,1)==BASIS%NUMBER_OF_XI) THEN
            IF(ANY(NUMBER_OF_GAUSS_XI<1)) CALL FLAG_ERROR("Invalid number of gauss values",ERR,ERROR,*999)
            BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI=NUMBER_OF_GAUSS_XI
            !Check the number of gauss points is sufficient for the interpolation order and flag a warning if not
            DO ni=1,BASIS%NUMBER_OF_XI
              SELECT CASE(BASIS%INTERPOLATION_XI(ni))
              CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<2) THEN
                  LOCAL_WARNING=TRIM(NUMBER_TO_VSTRING(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",ERR,ERROR))// &
                    & " Gauss points are insufficient for linear Lagrange interpolation"
                  CALL FLAG_WARNING(LOCAL_WARNING,ERR,ERROR,*999)
                ENDIF
              CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<2) THEN
                  LOCAL_WARNING=TRIM(NUMBER_TO_VSTRING(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",ERR,ERROR))//&
                    & " Gauss points are insufficient for quadratic Lagrange interpolation"
                  CALL FLAG_WARNING(LOCAL_WARNING,ERR,ERROR,*999)
                ENDIF
              CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<3) THEN
                  LOCAL_WARNING=TRIM(NUMBER_TO_VSTRING(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",ERR,ERROR))//&
                    & " Gauss points are insufficient for cubic Lagrange interpolation"
                  CALL FLAG_WARNING(LOCAL_WARNING,ERR,ERROR,*999)
                ENDIF
               CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<2) THEN
                  LOCAL_WARNING=TRIM(NUMBER_TO_VSTRING(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",ERR,ERROR))//&
                    & " Gauss points are insufficient for quadratic Hermite interpolation"
                  CALL FLAG_WARNING(LOCAL_WARNING,ERR,ERROR,*999)
                ENDIF
              CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<3) THEN
                  LOCAL_WARNING=TRIM(NUMBER_TO_VSTRING(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",ERR,ERROR))//&
                    & " Gauss points are insufficient for cubic Hermite interpolation"
                  CALL FLAG_WARNING(LOCAL_WARNING,ERR,ERROR,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="Interpolation xi value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_XI(ni),"*",ERR,ERROR))// &
                  & " is invalid for xi direction "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !xi
          ELSE
            LOCAL_ERROR="The size of the number of Gauss array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_GAUSS_XI,1),"*",ERR,ERROR))// &
              & ") does not match the number of xi directions ("// &
              & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//") for basis number "// &
              & TRIM(NUMBER_TO_VSTRING(BASIS%USER_NUMBER,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Quadrature basis is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET_PTR

  !
  !================================================================================================================================
  !
  
  !>Get the order of a quadrature for a basis quadrature identified by a pointer. \see OPENCMISS::CMISSBasisQuadratureOrderGet
  SUBROUTINE BASIS_QUADRATURE_ORDER_GET(BASIS,QUADRATURE_ORDER,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: QUADRATURE_ORDER !<On return, the quadrature order for the specified basis.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BASIS_QUADRATURE_ORDER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          QUADRATURE_ORDER=BASIS%QUADRATURE%GAUSS_ORDER
        ELSE
          CALL FLAG_ERROR("Quadrature basis is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Basis has not finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("BASIS_QUADRATURE_ORDER_GET")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_ORDER_GET",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_ORDER_GET")
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_ORDER_GET

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the order of a quadrature for a basis quadrature identified by a user number.
  SUBROUTINE BASIS_QUADRATURE_ORDER_SET_NUMBER(USER_NUMBER,ORDER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis
    INTEGER(INTG), INTENT(IN) :: ORDER !<The quadrature order to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_QUADRATURE_ORDER_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_ORDER_SET(BASIS,ORDER,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_QUADRATURE_ORDER_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_ORDER_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_ORDER_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_ORDER_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the order of a quadrature for a basis quadrature identified by a pointer. \see OPENCMISS::CMISSBasisQuadratureOrderSet
  SUBROUTINE BASIS_QUADRATURE_ORDER_SET_PTR(BASIS,ORDER,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: ORDER !<The quadrature order to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_QUADRATURE_ORDER_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          IF(BASIS%TYPE==BASIS_SIMPLEX_TYPE) THEN !Relax this i.e., use this to set gauss points in each direction for LHTP's???
            IF(ORDER>1.AND.ORDER<=5) THEN
              BASIS%QUADRATURE%GAUSS_ORDER=ORDER
            ELSE
              LOCAL_ERROR="An order value of "//TRIM(NUMBER_TO_VSTRING(ORDER,"*",ERR,ERROR))// &
                & " is invalid. You must specify and order between 1 and 5"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Can only set the quadrature order for simplex basis types",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Quadrature basis is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("BASIS_QUADRATURE_ORDER_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_ORDER_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_ORDER_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_ORDER_SET_PTR

  !
  !================================================================================================================================
  !
  
  !>get the quadrature type on a basis identified by a pointer. \see OPENCMISS::CMISSBasisQuadratureTypeGet
  SUBROUTINE BASIS_QUADRATURE_TYPE_GET(BASIS,QUADRATURE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: QUADRATURE_TYPE !<On return, the quadrature type to be get \see BASIS_ROUTINES_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BASIS_QUADRATURE_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          QUADRATURE_TYPE=BASIS%QUADRATURE%TYPE
        ELSE
          CALL FLAG_ERROR("Basis quadrature basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Basis has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_TYPE_GET")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_TYPE_GET",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_TYPE_GET")
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the quadrature type for a basis quadrature identified by a user number.
  SUBROUTINE BASIS_QUADRATURE_TYPE_SET_NUMBER(USER_NUMBER,TYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis
    INTEGER(INTG), INTENT(IN) :: TYPE !<The quadrature type to be set \see BASIS_ROUTINES_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_QUADRATURE_TYPE_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_TYPE_SET_PTR(BASIS,TYPE,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_QUADRATURE_TYPE_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_TYPE_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_TYPE_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the quadrature type on a basis identified by a pointer.
  SUBROUTINE BASIS_QUADRATURE_TYPE_SET_PTR(BASIS,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: TYPE !<The quadrature type to be set \see BASIS_ROUTINES_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_QUADRATURE_TYPE_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          SELECT CASE(TYPE)
          CASE(BASIS_GAUSS_LEGENDRE_QUADRATURE)
            BASIS%QUADRATURE%TYPE=BASIS_GAUSS_LEGENDRE_QUADRATURE
          CASE(BASIS_GAUSS_LAGUERRE_QUADRATURE)
            BASIS%QUADRATURE%TYPE=BASIS_GAUSS_LAGUERRE_QUADRATURE
            CALL FLAG_ERROR("Gauss Laguerre quadrature is not implemented",ERR,ERROR,*999)
          CASE(BASIS_GUASS_HERMITE_QUADRATURE)
            BASIS%QUADRATURE%TYPE=BASIS_GUASS_HERMITE_QUADRATURE
            CALL FLAG_ERROR("Gauss Hermite quadrature is not implemented",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Quadrature type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Basis quadrature basis is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_TYPE_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_TYPE_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_TYPE_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_TYPE_SET_PTR

  !
  !================================================================================================================================
  !

  !>Sets/changes the local face Gauss evaluation flag on a basis
  SUBROUTINE BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET(BASIS,FACE_GAUSS_EVALUATE,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    LOGICAL, INTENT(IN) :: FACE_GAUSS_EVALUATE !<face Gauss evaluation flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        BASIS%QUADRATURE%EVALUATE_FACE_GAUSS=FACE_GAUSS_EVALUATE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET")
    RETURN
999 CALL ERRORS("BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET",ERR,ERROR)
    CALL EXITS("BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET")
    RETURN 1

  END SUBROUTINE BASIS_QUADRATURE_LOCAL_FACE_GAUSS_EVALUATE_SET

  !
  !================================================================================================================================
  !
  
  !>Creates and initialises a simplex basis that has already been allocated BASIS_CREATE_START
  !>\see BASIS_ROUTINES::BASIS_CREATE_START
  SUBROUTINE BASIS_SIMPLEX_BASIS_CREATE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MAX_NUM_NODES,ni,nn,ns
    INTEGER(INTG), ALLOCATABLE :: NODES_IN_FACE(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_SIMPLEX_BASIS_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%TYPE==BASIS_SIMPLEX_TYPE) THEN
        BASIS%NUMBER_OF_XI_COORDINATES=BASIS%NUMBER_OF_XI+1 !Simplex bases have an additional area coordinate
        ALLOCATE(BASIS%INTERPOLATION_TYPE(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate INTERPOLATION_TYPE array",ERR,ERROR,*999)
        ALLOCATE(BASIS%INTERPOLATION_ORDER(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate INTERPOLATION_ORDER array",ERR,ERROR,*999)
        ALLOCATE(BASIS%NUMBER_OF_NODES_XIC(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NUMBER_OF_NODES_XIC array",ERR,ERROR,*999)
        BASIS%DEGENERATE=.FALSE.
        BASIS%NUMBER_OF_COLLAPSED_XI=0
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=3
          SELECT CASE(BASIS%INTERPOLATION_XI(1))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=2            
            BASIS%NUMBER_OF_NODES_XIC(2)=2            
            MAX_NUM_NODES=2
            BASIS%NUMBER_OF_NODES=2
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=3
            BASIS%NUMBER_OF_NODES_XIC(2)=3
            MAX_NUM_NODES=3
            BASIS%NUMBER_OF_NODES=3
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=4
            BASIS%NUMBER_OF_NODES_XIC(2)=4
            MAX_NUM_NODES=4
            BASIS%NUMBER_OF_NODES=4
          CASE DEFAULT 
            CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=6
          SELECT CASE(BASIS%INTERPOLATION_XI(2))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=2            
            BASIS%NUMBER_OF_NODES_XIC(2)=2            
            BASIS%NUMBER_OF_NODES_XIC(3)=2            
            MAX_NUM_NODES=2
            BASIS%NUMBER_OF_NODES=3
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=3
            BASIS%NUMBER_OF_NODES_XIC(2)=3
            BASIS%NUMBER_OF_NODES_XIC(3)=3
            MAX_NUM_NODES=3
            BASIS%NUMBER_OF_NODES=6
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=4
            BASIS%NUMBER_OF_NODES_XIC(2)=4
            BASIS%NUMBER_OF_NODES_XIC(3)=4
            MAX_NUM_NODES=4
            BASIS%NUMBER_OF_NODES=10
          CASE DEFAULT 
            CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
          END SELECT
        CASE(3)
          BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=11
          SELECT CASE(BASIS%INTERPOLATION_XI(3))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(4)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(4)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=2            
            BASIS%NUMBER_OF_NODES_XIC(2)=2            
            BASIS%NUMBER_OF_NODES_XIC(3)=2            
            BASIS%NUMBER_OF_NODES_XIC(4)=2            
            MAX_NUM_NODES=2
            BASIS%NUMBER_OF_NODES=4
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(4)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(4)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=3
            BASIS%NUMBER_OF_NODES_XIC(2)=3
            BASIS%NUMBER_OF_NODES_XIC(3)=3
            BASIS%NUMBER_OF_NODES_XIC(4)=3
            MAX_NUM_NODES=3
            BASIS%NUMBER_OF_NODES=10
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(4)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(4)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=4
            BASIS%NUMBER_OF_NODES_XIC(2)=4
            BASIS%NUMBER_OF_NODES_XIC(3)=4
            BASIS%NUMBER_OF_NODES_XIC(4)=4
            MAX_NUM_NODES=4
            BASIS%NUMBER_OF_NODES=20
          CASE DEFAULT 
            CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid number of xi directions",ERR,ERROR,*999)
        END SELECT
        
        ALLOCATE(BASIS%NODE_AT_COLLAPSE(BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node at collapse",ERR,ERROR,*999)
        BASIS%NODE_AT_COLLAPSE=.FALSE.
        
        ALLOCATE(BASIS%NODE_POSITION_INDEX(BASIS%NUMBER_OF_NODES,BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODE_POSITION_INDEX",ERR,ERROR,*999) 
        SELECT CASE(BASIS%NUMBER_OF_XI_COORDINATES)
        CASE(2)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,1,1),STAT=ERR)
        CASE(3)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES,1),STAT=ERR)
        CASE(4)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES),STAT=ERR)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
        END SELECT
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODE_POSITION_INDEX_INV",ERR,ERROR,*999)
        BASIS%NODE_POSITION_INDEX_INV=0
        
        !Determine the node position index and it's inverse
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=2
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,1)=2
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=3
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=2
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=3
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,1)=3
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=4
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX_INV(4,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=3
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX_INV(3,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=2
            BASIS%NODE_POSITION_INDEX(3,2)=3
            BASIS%NODE_POSITION_INDEX_INV(2,3,1,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=4
            BASIS%NODE_POSITION_INDEX_INV(1,4,1,1)=4
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid interpolation order",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=2
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,1)=3
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=3
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=3
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=3
            BASIS%NODE_POSITION_INDEX_INV(1,1,3,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=2
            BASIS%NODE_POSITION_INDEX(4,2)=2
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,1)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=1
            BASIS%NODE_POSITION_INDEX(5,2)=2
            BASIS%NODE_POSITION_INDEX(5,3)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,2,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=2
            BASIS%NODE_POSITION_INDEX(6,2)=1
            BASIS%NODE_POSITION_INDEX(6,3)=2
            BASIS%NODE_POSITION_INDEX_INV(2,1,2,1)=6
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=4
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX_INV(4,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=4
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX_INV(1,4,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=4
            BASIS%NODE_POSITION_INDEX_INV(1,1,4,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=3
            BASIS%NODE_POSITION_INDEX(4,2)=2
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX_INV(3,2,1,1)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=2
            BASIS%NODE_POSITION_INDEX(5,2)=3
            BASIS%NODE_POSITION_INDEX(5,3)=1
            BASIS%NODE_POSITION_INDEX_INV(2,3,1,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=1
            BASIS%NODE_POSITION_INDEX(6,2)=3
            BASIS%NODE_POSITION_INDEX(6,3)=2
            BASIS%NODE_POSITION_INDEX_INV(1,3,2,1)=6
            !Node 7
            BASIS%NODE_POSITION_INDEX(7,1)=1
            BASIS%NODE_POSITION_INDEX(7,2)=2
            BASIS%NODE_POSITION_INDEX(7,3)=3
            BASIS%NODE_POSITION_INDEX_INV(1,2,3,1)=7
            !Node 8
            BASIS%NODE_POSITION_INDEX(8,1)=2
            BASIS%NODE_POSITION_INDEX(8,2)=1
            BASIS%NODE_POSITION_INDEX(8,3)=3
            BASIS%NODE_POSITION_INDEX_INV(2,1,3,1)=8
            !Node 9
            BASIS%NODE_POSITION_INDEX(9,1)=3
            BASIS%NODE_POSITION_INDEX(9,2)=1
            BASIS%NODE_POSITION_INDEX(9,3)=2
            BASIS%NODE_POSITION_INDEX_INV(3,1,2,1)=9
            !Node 10
            BASIS%NODE_POSITION_INDEX(10,1)=2
            BASIS%NODE_POSITION_INDEX(10,2)=2
            BASIS%NODE_POSITION_INDEX(10,3)=2
            BASIS%NODE_POSITION_INDEX_INV(2,2,2,1)=10
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid interpolation order",ERR,ERROR,*999)
          END SELECT
        CASE(3)
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=2
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX(1,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX(2,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=2
            BASIS%NODE_POSITION_INDEX(3,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=1
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX(4,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,1,2)=4

            ALLOCATE(NODES_IN_FACE(12),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES_IN_FACE",ERR,ERROR,*999) 
            NODES_IN_FACE(:)=(/2,3,4,1,3,4,1,2,4,1,2,3/) !12 Nodes

          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=3
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX(1,4)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=3
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX(2,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=3
            BASIS%NODE_POSITION_INDEX(3,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,1,3,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=1
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX(4,4)=3
            BASIS%NODE_POSITION_INDEX_INV(1,1,1,3)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=2
            BASIS%NODE_POSITION_INDEX(5,2)=2
            BASIS%NODE_POSITION_INDEX(5,3)=1
            BASIS%NODE_POSITION_INDEX(5,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=2
            BASIS%NODE_POSITION_INDEX(6,2)=1
            BASIS%NODE_POSITION_INDEX(6,3)=2
            BASIS%NODE_POSITION_INDEX(6,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,2,1)=6
            !Node 7
            BASIS%NODE_POSITION_INDEX(7,1)=2
            BASIS%NODE_POSITION_INDEX(7,2)=1
            BASIS%NODE_POSITION_INDEX(7,3)=1
            BASIS%NODE_POSITION_INDEX(7,4)=2
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,2)=7
            !Node 8
            BASIS%NODE_POSITION_INDEX(8,1)=1
            BASIS%NODE_POSITION_INDEX(8,2)=2
            BASIS%NODE_POSITION_INDEX(8,3)=2
            BASIS%NODE_POSITION_INDEX(8,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,2,1)=8
            !Node 9
            BASIS%NODE_POSITION_INDEX(9,1)=1
            BASIS%NODE_POSITION_INDEX(9,2)=1
            BASIS%NODE_POSITION_INDEX(9,3)=2
            BASIS%NODE_POSITION_INDEX(9,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,2)=9
            !Node 10
            BASIS%NODE_POSITION_INDEX(10,1)=1
            BASIS%NODE_POSITION_INDEX(10,2)=2
            BASIS%NODE_POSITION_INDEX(10,3)=1
            BASIS%NODE_POSITION_INDEX(10,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,2)=10

            ALLOCATE(NODES_IN_FACE(24),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES_IN_FACE",ERR,ERROR,*999) 
            NODES_IN_FACE(:)=(/2,3,4,8,9,10,1,3,4,6,9,7,&
                               &1,2,4,5,10,7,1,2,3,5,8,6/) !24 Nodes

          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=4
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX(1,4)=1
            BASIS%NODE_POSITION_INDEX_INV(4,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=4
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX(2,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,4,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=4
            BASIS%NODE_POSITION_INDEX(3,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,1,4,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=1
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX(4,4)=4
            BASIS%NODE_POSITION_INDEX_INV(1,1,1,4)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=3
            BASIS%NODE_POSITION_INDEX(5,2)=2
            BASIS%NODE_POSITION_INDEX(5,3)=1
            BASIS%NODE_POSITION_INDEX(5,4)=1
            BASIS%NODE_POSITION_INDEX_INV(3,2,1,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=2
            BASIS%NODE_POSITION_INDEX(6,2)=3
            BASIS%NODE_POSITION_INDEX(6,3)=1
            BASIS%NODE_POSITION_INDEX(6,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,3,1,1)=6
            !Node 7
            BASIS%NODE_POSITION_INDEX(7,1)=3
            BASIS%NODE_POSITION_INDEX(7,2)=1
            BASIS%NODE_POSITION_INDEX(7,3)=2
            BASIS%NODE_POSITION_INDEX(7,4)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,2,1)=7
            !Node 8
            BASIS%NODE_POSITION_INDEX(8,1)=2
            BASIS%NODE_POSITION_INDEX(8,2)=1
            BASIS%NODE_POSITION_INDEX(8,3)=3
            BASIS%NODE_POSITION_INDEX(8,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,3,1)=8
            !Node 9
            BASIS%NODE_POSITION_INDEX(9,1)=3
            BASIS%NODE_POSITION_INDEX(9,2)=1
            BASIS%NODE_POSITION_INDEX(9,3)=1
            BASIS%NODE_POSITION_INDEX(9,4)=2
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,2)=9
            !Node 10
            BASIS%NODE_POSITION_INDEX(10,1)=2
            BASIS%NODE_POSITION_INDEX(10,2)=1
            BASIS%NODE_POSITION_INDEX(10,3)=1
            BASIS%NODE_POSITION_INDEX(10,4)=3
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,3)=10
            !Node 11
            BASIS%NODE_POSITION_INDEX(11,1)=1
            BASIS%NODE_POSITION_INDEX(11,2)=3
            BASIS%NODE_POSITION_INDEX(11,3)=2
            BASIS%NODE_POSITION_INDEX(11,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,3,2,1)=11
            !Node 12
            BASIS%NODE_POSITION_INDEX(12,1)=1
            BASIS%NODE_POSITION_INDEX(12,2)=2
            BASIS%NODE_POSITION_INDEX(12,3)=3
            BASIS%NODE_POSITION_INDEX(12,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,3,1)=12
            !Node 13
            BASIS%NODE_POSITION_INDEX(13,1)=1
            BASIS%NODE_POSITION_INDEX(13,2)=1
            BASIS%NODE_POSITION_INDEX(13,3)=3
            BASIS%NODE_POSITION_INDEX(13,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,3,2)=13
            !Node 14
            BASIS%NODE_POSITION_INDEX(14,1)=1
            BASIS%NODE_POSITION_INDEX(14,2)=1
            BASIS%NODE_POSITION_INDEX(14,3)=2
            BASIS%NODE_POSITION_INDEX(14,4)=3
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,3)=14
            !Node 15
            BASIS%NODE_POSITION_INDEX(15,1)=1
            BASIS%NODE_POSITION_INDEX(15,2)=3
            BASIS%NODE_POSITION_INDEX(15,3)=1
            BASIS%NODE_POSITION_INDEX(15,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,2)=15
            !Node 16
            BASIS%NODE_POSITION_INDEX(16,1)=1
            BASIS%NODE_POSITION_INDEX(16,2)=2
            BASIS%NODE_POSITION_INDEX(16,3)=1
            BASIS%NODE_POSITION_INDEX(16,4)=3
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,3)=16
            !Node 17
            BASIS%NODE_POSITION_INDEX(17,1)=2
            BASIS%NODE_POSITION_INDEX(17,2)=2
            BASIS%NODE_POSITION_INDEX(17,3)=2
            BASIS%NODE_POSITION_INDEX(17,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,2,2,1)=17
            !Node 18
            BASIS%NODE_POSITION_INDEX(18,1)=2
            BASIS%NODE_POSITION_INDEX(18,2)=2
            BASIS%NODE_POSITION_INDEX(18,3)=1
            BASIS%NODE_POSITION_INDEX(18,4)=2
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,2)=18
            !Node 19
            BASIS%NODE_POSITION_INDEX(19,1)=2
            BASIS%NODE_POSITION_INDEX(19,2)=1
            BASIS%NODE_POSITION_INDEX(19,3)=2
            BASIS%NODE_POSITION_INDEX(19,4)=2
            BASIS%NODE_POSITION_INDEX_INV(2,1,2,2)=19
            !Node 20
            BASIS%NODE_POSITION_INDEX(20,1)=1
            BASIS%NODE_POSITION_INDEX(20,2)=2
            BASIS%NODE_POSITION_INDEX(20,3)=2
            BASIS%NODE_POSITION_INDEX(20,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,2,2)=20

            ALLOCATE(NODES_IN_FACE(40),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES_IN_FACE",ERR,ERROR,*999) 
            NODES_IN_FACE(:)=(/2,3,4,11,12,13,14,16,15,20,1,3,4,7,8,13,14,10,9,&
                               &19,1,2,4,5,6,15,16,10,9,18,1,2,3,5,6,14,12,8,7,17/) !40 nodes

          CASE DEFAULT
            CALL FLAG_ERROR("Invalid interpolation order",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid number of xi directions",ERR,ERROR,*999)
        END SELECT
        !Calculate the maximum number of derivatives (1 for simplex bases) and the number of element parameters
        BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES=1
        BASIS%NUMBER_OF_ELEMENT_PARAMETERS=BASIS%NUMBER_OF_NODES
        !Now set up the number of derivatives and derivative order index
        ALLOCATE(BASIS%NUMBER_OF_DERIVATIVES(BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NUMBER_OF_DERIVATIVES",ERR,ERROR,*999)
        ALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES,BASIS%NUMBER_OF_XI), &
          & STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DERIVATIVE_ORDER_INDEX",ERR,ERROR,*999)
        ALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX_INV(FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,BASIS%NUMBER_OF_NODES), &
          & STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DERIVATIVE_ORDER_INDEX_INV",ERR,ERROR,*999)
        ALLOCATE(BASIS%PARTIAL_DERIVATIVE_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate PARTIAL_DERIVATIVE_INDEX",ERR,ERROR,*999)
        ALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT_PARAMETER_INDEX",ERR,ERROR,*999)
        ALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX_INV(2,BASIS%NUMBER_OF_ELEMENT_PARAMETERS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT_PARAMETER_INDEX_INV",ERR,ERROR,*999)
        !Set the derivative order index and its inverse, the element parameter index and the partial derivative index.
        ns=0
        BASIS%DERIVATIVE_ORDER_INDEX_INV=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          BASIS%NUMBER_OF_DERIVATIVES(nn)=1
          DO ni=1,BASIS%NUMBER_OF_XI
            BASIS%DERIVATIVE_ORDER_INDEX(1,nn,ni)=1
          ENDDO !ni
          ns=ns+1
          BASIS%ELEMENT_PARAMETER_INDEX(1,nn)=ns
          BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ns)=nn
          BASIS%ELEMENT_PARAMETER_INDEX_INV(2,ns)=1
          BASIS%PARTIAL_DERIVATIVE_INDEX(1,nn)=NO_PART_DERIV
          BASIS%DERIVATIVE_ORDER_INDEX_INV(BASIS%DERIVATIVE_ORDER_INDEX(1,nn,1),1,1,nn)=1
        ENDDO !nn
      
        !Set up the line and face information
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          BASIS%NUMBER_OF_LOCAL_LINES=1
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local line",ERR,ERROR,*999)
          ALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line xi direction",ERR,ERROR,*999)
          BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(BASIS%NUMBER_OF_NODES_XIC(1),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local line",ERR,ERROR,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(BASIS%NUMBER_OF_NODES_XIC(1),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local line",ERR,ERROR,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          !Set the line values
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=3
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,1)=4
          CASE DEFAULT 
            LOCAL_ERROR="Interpolation order "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(1),"*",ERR,ERROR))// &
              & " is invalid for a simplex basis type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(2)
          !Allocate and calculate the lines
          !Simplex hence three local lines
          BASIS%NUMBER_OF_LOCAL_LINES=3
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local line.",ERR,ERROR,*999)
          ALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line xi direction",ERR,ERROR,*999)
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local line.",ERR,ERROR,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local line.",ERR,ERROR,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(BASIS%LOCAL_XI_NORMAL(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line normal.",ERR,ERROR,*999)
          !Set the line values
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
            BASIS%LOCAL_XI_NORMAL(1)=3
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(2)=2
            BASIS%LOCAL_XI_NORMAL(2)=2
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(3)=3
            BASIS%LOCAL_XI_NORMAL(3)=1
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=2
            BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
            BASIS%LOCAL_XI_NORMAL(1)=3
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(2)=2
            BASIS%LOCAL_XI_NORMAL(2)=2
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(3)=3
            BASIS%LOCAL_XI_NORMAL(3)=1
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,1)=2
            BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
            BASIS%LOCAL_XI_NORMAL(1)=3
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,2)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(2)=2
            BASIS%LOCAL_XI_NORMAL(2)=2
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,3)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(3)=3
            BASIS%LOCAL_XI_NORMAL(3)=1
          CASE DEFAULT 
            LOCAL_ERROR="Interpolation order "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(1),"*",ERR,ERROR))// &
              & " is invalid for a simplex basis type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(3)
          BASIS%NUMBER_OF_LOCAL_LINES=6
          BASIS%NUMBER_OF_LOCAL_FACES=4

          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local line.",ERR,ERROR,*999)
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of nodes in local face.",ERR,ERROR,*999)

          ALLOCATE(BASIS%LOCAL_LINE_XI_DIRECTION(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line xi direction.",ERR,ERROR,*999)
          ALLOCATE(BASIS%LOCAL_FACE_XI_DIRECTION(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local face xi direction.",ERR,ERROR,*999)

          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local line.",ERR,ERROR,*999)
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node numbers in local face.",ERR,ERROR,*999)

          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local line.",ERR,ERROR,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(10,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivative numbers in local face.",ERR,ERROR,*999)

          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV

          ALLOCATE(BASIS%LOCAL_XI_NORMAL(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local line normal.",ERR,ERROR,*999)

          !Set the line values
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(2)=1
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(3)=1
            !Line 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,4)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(4)=2
            !Line 5
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,5)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(5)=2
            !Line 6
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(6)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,6)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(6)=3
            !Face 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,1)=4
            BASIS%LOCAL_FACE_XI_DIRECTION(1)=1
            BASIS%LOCAL_XI_NORMAL(1)=1
            !Face 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,2)=4
            BASIS%LOCAL_FACE_XI_DIRECTION(2)=2
            BASIS%LOCAL_XI_NORMAL(2)=2
            !Face 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(3)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,3)=4
            BASIS%LOCAL_FACE_XI_DIRECTION(3)=3
            BASIS%LOCAL_XI_NORMAL(3)=3
            !Face 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,4)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,4)=3
            BASIS%LOCAL_FACE_XI_DIRECTION(4)=4
            BASIS%LOCAL_XI_NORMAL(4)=4
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=2
            BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(2)=1
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(3)=1
            !Line 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,4)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,4)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(4)=2
            !Line 5
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(5)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,5)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,5)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(5)=2
            !Line 6
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,6)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,6)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(6)=3
            !Face 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(1)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,1)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,1)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,1)=10
            BASIS%LOCAL_FACE_XI_DIRECTION(1)=1
            BASIS%LOCAL_XI_NORMAL(1)=1
            !Face 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,2)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,2)=7
            BASIS%LOCAL_FACE_XI_DIRECTION(2)=2
            BASIS%LOCAL_XI_NORMAL(2)=2
            !Face 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(3)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,3)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,3)=7
            BASIS%LOCAL_FACE_XI_DIRECTION(3)=3
            BASIS%LOCAL_XI_NORMAL(3)=3
            !Face 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(4)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,4)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,4)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,4)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,4)=6
            BASIS%LOCAL_FACE_XI_DIRECTION(4)=4
            BASIS%LOCAL_XI_NORMAL(4)=4
           CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,1)=2
            BASIS%LOCAL_LINE_XI_DIRECTION(1)=1
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,2)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(2)=1
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,3)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(3)=1
            !Line 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(4)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,4)=11
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,4)=12
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,4)=3
            BASIS%LOCAL_LINE_XI_DIRECTION(4)=2
            !Line 5
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(5)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,5)=15
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,5)=16
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,5)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(5)=2
            !Line 6
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(6)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,6)=13
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,6)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,6)=4
            BASIS%LOCAL_LINE_XI_DIRECTION(6)=3
            !Face 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(1)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,1)=11
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,1)=12
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,1)=13
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,1)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,1)=16
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,1)=15
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,1)=20
            BASIS%LOCAL_FACE_XI_DIRECTION(1)=1
            BASIS%LOCAL_XI_NORMAL(1)=1
             !Face 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(2)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,2)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,2)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,2)=13
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,2)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,2)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,2)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,2)=19
            BASIS%LOCAL_FACE_XI_DIRECTION(2)=2
            BASIS%LOCAL_XI_NORMAL(2)=2
            !Face 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,3)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,3)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,3)=15
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,3)=16
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,3)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,3)=18
            BASIS%LOCAL_FACE_XI_DIRECTION(3)=3
            BASIS%LOCAL_XI_NORMAL(3)=3
            !Face 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(4)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,4)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,4)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,4)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,4)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,4)=12
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,4)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,4)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,4)=17
            BASIS%LOCAL_FACE_XI_DIRECTION(4)=4
            BASIS%LOCAL_XI_NORMAL(4)=4
           CASE DEFAULT
            LOCAL_ERROR="Interpolation order "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(1),"*",ERR,ERROR))// &
              & " is invalid for a simplex basis type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid number of xi directions.",ERR,ERROR,*999)
        END SELECT
        
        CALL BASIS_QUADRATURE_CREATE(BASIS,ERR,ERROR,*999)
        
      ELSE
        CALL FLAG_ERROR("Basis is not a simplex basis.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BASIS_SIMPLEX_BASIS_CREATE")
    RETURN
999 CALL ERRORS("BASIS_SIMPLEX_BASIS_CREATE",ERR,ERROR)
    CALL EXITS("BASIS_SIMPLEX_BASIS_CREATE")
    RETURN 1
  END SUBROUTINE BASIS_SIMPLEX_BASIS_CREATE

  !
  !================================================================================================================================
  !

  !>Creates and initialises a simplex basis family that has already been allocated by BASIS_CREATE_START.
  !> \see BASIS_ROUTINES::BASIS_SIMPLEX_BASIS_CREATE
  SUBROUTINE BASIS_SIMPLEX_FAMILY_CREATE(BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,ni,ni2,FACE_XI(2),FACE_XI2(2)
    LOGICAL :: LINE_BASIS_DONE,FACE_BASIS_DONE
    TYPE(BASIS_TYPE), POINTER :: NEW_SUB_BASIS
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_SUB_BASIS)

    CALL ENTERS("BASIS_SIMPLEX_FAMILY_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      !Create the main (parent) basis
      CALL BASIS_SIMPLEX_BASIS_CREATE(BASIS,ERR,ERROR,*999)
      IF(BASIS%NUMBER_OF_XI>1) THEN
        !Create the line bases as sub-basis types
!        ALLOCATE(BASIS%LINE_BASES(BASIS%NUMBER_OF_XI),STAT=ERR)
        IF (BASIS%NUMBER_OF_XI .EQ. 2) THEN
          ALLOCATE(BASIS%LINE_BASES(BASIS%NUMBER_OF_XI+1),STAT=ERR)
        ELSE
          ALLOCATE(BASIS%LINE_BASES(BASIS%NUMBER_OF_XI),STAT=ERR)
        ENDIF
 
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis line bases",ERR,ERROR,*999)
        DO ni=1,BASIS%NUMBER_OF_XI
          LINE_BASIS_DONE=.FALSE.
          NULLIFY(NEW_SUB_BASIS)
          DO ni2=1,ni-1
            IF(BASIS%INTERPOLATION_XI(ni2)==BASIS%INTERPOLATION_XI(ni).AND. &
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni2)==BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)) THEN
              LINE_BASIS_DONE=.TRUE.
              EXIT
            ENDIF
          ENDDO !ni2
          IF(LINE_BASIS_DONE) THEN
            BASIS%LINE_BASES(ni)%PTR=>BASIS%LINE_BASES(ni2)%PTR
          ELSE
            !Create the new sub-basis
            CALL BASIS_SUB_BASIS_CREATE(BASIS,1,(/ni/),NEW_SUB_BASIS,ERR,ERROR,*999)
            !Fill in the basis information
            CALL BASIS_SIMPLEX_BASIS_CREATE(NEW_SUB_BASIS,ERR,ERROR,*999)
            BASIS%LINE_BASES(ni)%PTR=>NEW_SUB_BASIS
          ENDIF
        ENDDO !ni

        IF (BASIS%NUMBER_OF_XI .EQ. 2) THEN
          BASIS%LINE_BASES(BASIS%NUMBER_OF_XI+1)%PTR=>BASIS%LINE_BASES(BASIS%NUMBER_OF_XI)%PTR
        ENDIF
        
        IF(BASIS%NUMBER_OF_XI>2) THEN
          !Set up face basis functions
          ALLOCATE(BASIS%FACE_BASES(BASIS%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis face bases",ERR,ERROR,*999)
          DO ni=1,BASIS%NUMBER_OF_XI
            !Determine the face xi directions that lie in this xi direction
            FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
            FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
            FACE_BASIS_DONE=.FALSE.
            NULLIFY(NEW_SUB_BASIS)
            DO ni2=1,ni-1
              FACE_XI2(1)=OTHER_XI_DIRECTIONS3(ni2,2,1)
              FACE_XI2(2)=OTHER_XI_DIRECTIONS3(ni2,3,1)
              IF(BASIS%INTERPOLATION_XI(FACE_XI2(1))==BASIS%INTERPOLATION_XI(FACE_XI(1)).AND. &
                & BASIS%INTERPOLATION_XI(FACE_XI2(2))==BASIS%INTERPOLATION_XI(FACE_XI(2)).AND. &
                & BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI2(1))==BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(1)).AND. &
                & BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI2(2))==BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(1))) THEN
                FACE_BASIS_DONE=.TRUE.
                EXIT
              ENDIF
            ENDDO !ni2
            IF(FACE_BASIS_DONE) THEN
              BASIS%FACE_BASES(ni)%PTR=>BASIS%FACE_BASES(ni2)%PTR
            ELSE
              !Create the new sub-basis
              CALL BASIS_SUB_BASIS_CREATE(BASIS,2,(/FACE_XI(1),FACE_XI(2)/),NEW_SUB_BASIS,ERR,ERROR,*999)
              !Fill in the basis information
              CALL BASIS_SIMPLEX_BASIS_CREATE(NEW_SUB_BASIS,ERR,ERROR,*999)
              NEW_SUB_BASIS%LINE_BASES(1)%PTR=>BASIS%LINE_BASES(FACE_XI(1))%PTR
              NEW_SUB_BASIS%LINE_BASES(2)%PTR=>BASIS%LINE_BASES(FACE_XI(2))%PTR
              BASIS%FACE_BASES(ni)%PTR=>NEW_SUB_BASIS
            ENDIF            
          ENDDO !ni
        ELSE
          ALLOCATE(BASIS%FACE_BASES(1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis face bases",ERR,ERROR,*999)
          BASIS%FACE_BASES(1)%PTR=>BASIS
        ENDIF
      ELSE
        ALLOCATE(BASIS%LINE_BASES(1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis line bases",ERR,ERROR,*999)
        BASIS%LINE_BASES(1)%PTR=>BASIS
        NULLIFY(BASIS%FACE_BASES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_SIMPLEX_FAMILY_CREATE")
    RETURN
999 IF(ASSOCIATED(NEW_SUB_BASIS)) CALL BASIS_FAMILY_DESTROY(NEW_SUB_BASIS%USER_NUMBER,NEW_SUB_BASIS%FAMILY_NUMBER, &
      & DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("BASIS_SIMPLEX_FAMILY_CREATE",ERR,ERROR)
    CALL EXITS("BASIS_SIMPLEX_FAMILY_CREATE")
    RETURN 1
  END SUBROUTINE BASIS_SIMPLEX_FAMILY_CREATE

  !
  !================================================================================================================================
  !

  !>Evaluates a simplex basis function and its derivatives with respect to external \f$\mathbf{\xi}\f$ coordinates.
  !>For Simplex line elements there are two area coordinates which are a function of \f$\xi_1\f$ : \f$L_1 = 1 - \xi_1\f$ and
  !>\f$L_2 = \xi_1 - 1\f$.The derivatives wrt to external coordinates are then given by \f$\frac{\partial\mathbf{N}}{\partial\xi_1}=
  !>\frac{\partial\mathbf(x)}{\partial L_2}-\frac{\partial \mathbf{N}}{\partial L_1}\f$ and \f$\frac{\partial^2\mathbf{N}}{
  !>\partial \xi_1^2} = \frac{\partial^2\mathbf{N}}{\partial L_1^2}-2\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_2}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_2^2}\f$.
  !>For Simplex triangle elements there are three area coordinates which are a function of \f$\xi_1\f$ and
  !>\f$\xi_2\f$ : \f$L_1 = 1 - \xi_1\f$, \f$L_2 = 1 - \xi_2\f$ and \f$L_3=\xi_1 + \xi_2 - 1\f$. The derivatives wrt to external
  !>coordinates are then given by \f$\frac{\partial \mathbf{N}}{\partial\xi_1}=\frac{\partial\mathbf(N)}{\partial L_3}-
  !>\frac{\partial \mathbf{N}}{\partial L_1}\f$, \f$\frac{\partial \mathbf{N}}{\partial\xi_2}=\frac{\partial\mathbf(x)}{
  !>\partial L_3}-\frac{\partial \mathbf{N}}{\partial L_2}\f$, \f$\frac{\partial^2\mathbf{N}}{\partial \xi_1^2} =
  !>\frac{\partial^2\mathbf{N}}{\partial L_1^2}-2\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_3}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_3^2}\f$, \f$\frac{\partial^2\mathbf{N}}{\partial \xi_2^2} =
  !>\frac{\partial^2\mathbf{N}}{\partial L_2^2}-2\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_3}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_3^2}\f$ and \f$\frac{\partial^2\mathbf{N}}{\partial \xi_1 \partial \xi_2} =
  !>\frac{\partial^2\mathbf{N}}{\partial L_3^2}-\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_3}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_3}+\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_2}\f$.
  !>For Simplex tetrahedral elements there are four area coordinates which are a function of \f$\xi_1\f$,\f$\xi_2\f$ and
  !>\f$\xi_3\f$ : \f$L_1 = 1 - \xi_1\f$, \f$L_2 = 1 - \xi_2\f$, \f$L_3 = 1 - \xi_3\f$ and
  !>\f$L_4 = \xi_1 + \xi_2 + \xi_3 - 1\f$. The derivatives wrt to external coordinates are then given by
  !>\f$\frac{\partial \mathbf{N}}{\partial\xi_1}=\frac{\partial\mathbf(x)}{\partial L_4}-
  !>\frac{\partial \mathbf{N}}{\partial L_1}\f$,
  !>\f$\frac{\partial \mathbf{N}}{\partial\xi_2}=\frac{\partial\mathbf(x)}{\partial L_4}-
  !>\frac{\partial \mathbf{N}}{\partial L_2}\f$,
  !>\f$\frac{\partial \mathbf{N}}{\partial\xi_3}=\frac{\partial\mathbf(x)}{\partial L_4}-
  !>\frac{\partial \mathbf{N}}{\partial L_3}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial \xi_1^2} = \frac{\partial^2\mathbf{N}}{\partial L_1^2}-
  !>2\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_4}+\frac{\partial^2\mathbf{N}}{\partial L_4^2}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial \xi_2^2} = \frac{\partial^2\mathbf{N}}{\partial L_2^2}-
  !>2\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_4}+\frac{\partial^2\mathbf{N}}{\partial L_4^2}\f$
  !>\f$\frac{\partial^2\mathbf{N}}{\partial \xi_3^2} = \frac{\partial^2\mathbf{N}}{\partial L_3^2}-
  !>2\frac{\partial^2\mathbf{N}}{\partial L_3 \partial L_4}+\frac{\partial^2\mathbf{N}}{\partial L_4^2}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial\xi_1\partial \xi_2}=\frac{\partial^2\mathbf{N}}{\partial L_4^2}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_4}-\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_4}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_2}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial\xi_1\partial\xi_3}=\frac{\partial^2\mathbf{N}}{\partial L_4^2}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_4}-\frac{\partial^2\mathbf{N}}{\partial L_3 \partial L_4}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_3}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial\xi_2\partial\xi_3}=\frac{\partial^2\mathbf{N}}{\partial L_4^2}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_4}-\frac{\partial^2\mathbf{N}}{\partial L_3 \partial L_4}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_3}\f$ and
  !>\f$\frac{\partial^3\mathbf{N}}{\partial \xi_1 \partial \xi_2 \partial \xi_3} = \frac{\partial^3\mathbf{N}}{\partial L_4^3}-
  !>\frac{\partial^3\mathbf{N}}{\partial L_1 \partial L_4^2}-\frac{\partial^3\mathbf{N}}{\partial L_2 \partial L_4^2}-
  !>\frac{\partial^3\mathbf{N}}{\partial L_3 \partial L_4^2}+\frac{\partial^3\mathbf{N}}{\partial L_1 \partial 2 \partial L_4}+
  !>\frac{\partial^3\mathbf{N}}{\partial L_1 \partial L_3 \partial L_4}+\frac{\partial^3\mathbf{N}}{\partial L_2 \partial L_3
  !>\partial L_4}-\frac{\partial^3\mathbf{N}}{\partial L_1 \partial L_2 \partial L_3}\f$.
  FUNCTION BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,NODE_NUMBER,PARTIAL_DERIV_INDEX,XL,ERR,ERROR)
    
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The node number defines the actual basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index in Xi space of the basis to evaluate.
    REAL(DP), INTENT(IN) :: XL(:) !<XL(nic). The area coordinates to evaluate the basis function at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_SIMPLEX_BASIS_EVALUATE !<On return the evaluated basis function
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_SIMPLEX_BASIS_EVALUATE",ERR,ERROR,*999)
    
    BASIS_SIMPLEX_BASIS_EVALUATE=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%TYPE==BASIS_SIMPLEX_TYPE) THEN
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          SELECT CASE(PARTIAL_DERIV_INDEX)
          CASE(NO_PART_DERIV)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,NO_PART_DERIV,XL,ERR,ERROR)
          CASE(PART_DERIV_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S1,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S2,XL,ERR,ERROR)
          CASE DEFAULT
            LOCAL_ERROR="The specified partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid for a Simplex line basis."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(PARTIAL_DERIV_INDEX)
          CASE(NO_PART_DERIV)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,NO_PART_DERIV,XL,ERR,ERROR)
          CASE(PART_DERIV_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S1,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,ERR,ERROR)
          CASE(PART_DERIV_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2,XL,ERR,ERROR)
          CASE(PART_DERIV_S2_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S2,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2,XL,ERR,ERROR)
          CASE DEFAULT
            LOCAL_ERROR="The specified partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid for a Simplex triangle basis."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(3)
          SELECT CASE(PARTIAL_DERIV_INDEX)
          CASE(NO_PART_DERIV)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,NO_PART_DERIV,XL,ERR,ERROR)
          CASE(PART_DERIV_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S1)
             BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S1,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,ERR,ERROR)
          CASE(PART_DERIV_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2,XL,ERR,ERROR)
          CASE(PART_DERIV_S2_S2)
             BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S2,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2,XL,ERR,ERROR)
          CASE(PART_DERIV_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3,XL,ERR,ERROR)
          CASE(PART_DERIV_S3_S3)
             BASIS_SIMPLEX_BASIS_EVALUATE= &
               BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3,XL,ERR,ERROR)
          CASE(PART_DERIV_S2_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3,XL,ERR,ERROR)
          CASE(PART_DERIV_S1_S2_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3_S4,XL,ERR,ERROR)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2_S3,XL,ERR,ERROR)
          CASE DEFAULT
            LOCAL_ERROR="The specified partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid for a Simplex tetrahedra basis."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="Invalid number of Xi coordinates. The number of xi coordinates for this basis is "// &
            & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" which should be between 1 and 3."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(ERR/=0) GOTO 999
      ELSE
        CALL FLAG_ERROR("Basis is not a simplex basis.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_SIMPLEX_BASIS_EVALUATE.")
    RETURN
999 CALL ERRORS("BASIS_SIMPLEX_BASIS_EVALUATE.",ERR,ERROR)
    CALL EXITS("BASIS_SIMPLEX_BASIS_EVALUATE.")
    RETURN 
  END FUNCTION BASIS_SIMPLEX_BASIS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates partial derivatives of a simplex basis function with respect to area coordinates.
  FUNCTION BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PARTIAL_DERIV_INDEX,XL,ERR,ERROR)
    
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The node number defines the actual basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index in area coordinates of the basis to evaluate.
    REAL(DP), INTENT(IN) :: XL(:) !<XL(nic). The area coordinates to evaluate the basis function at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE !<On return the evaluated basis function
    !Local variables
    INTEGER(INTG) :: nic
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE",ERR,ERROR,*999)
    
    BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=1.0_DP
    IF(ASSOCIATED(BASIS)) THEN      
      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES       
        SELECT CASE(BASIS%INTERPOLATION_ORDER(nic))
        CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
          BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE* &
            & SIMPLEX_LINEAR_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,nic), &
            & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,nic),XL(nic),ERR,ERROR)
        CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
          BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE* &
            & SIMPLEX_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,nic), &
            & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,nic),XL(nic),ERR,ERROR)
        CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
          BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE* &
            & SIMPLEX_CUBIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,nic), &
            & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,nic),XL(nic),ERR,ERROR)
        CASE DEFAULT
          LOCAL_ERROR="Interpolation order value "//TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(nic),"*",ERR,ERROR))// &
            & " for xi coordinate direction "//TRIM(NUMBER_TO_VSTRING(nic,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(ERR/=0) GOTO 999
      ENDDO !nic
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE")
    RETURN
999 CALL ERRORS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE",ERR,ERROR)
    CALL EXITS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE")
    RETURN 
  END FUNCTION BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE

  !
  !================================================================================================================================
  !

  !>Creates a sub-basis on a parent basis.
  SUBROUTINE BASIS_SUB_BASIS_CREATE(PARENT_BASIS,NUMBER_OF_XI,XI_DIRECTIONS,SUB_BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: PARENT_BASIS !<A pointer to the parent basis
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI !<The number of Xi directions to create
    INTEGER(INTG), INTENT(IN) :: XI_DIRECTIONS(:) !<XI_DIRECTIONS(ni). Gives the Xi direction indices of the parent basis which are used to create the sub-basis
    TYPE(BASIS_TYPE), POINTER :: SUB_BASIS !<On return, a pointer to the created sub-basis. The pointer must be NULL on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: basis_idx,ni,NUMBER_COLLAPSED,NUMBER_END_COLLAPSED
    TYPE(BASIS_TYPE), POINTER :: NEW_SUB_BASIS    
    TYPE(BASIS_PTR_TYPE), POINTER :: NEW_SUB_BASES(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    NULLIFY(NEW_SUB_BASIS)
    NULLIFY(NEW_SUB_BASES)
    
    CALL ENTERS("BASIS_SUB_BASIS_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(PARENT_BASIS)) THEN
      IF(ASSOCIATED(SUB_BASIS)) THEN
        CALL FLAG_ERROR("The sub-basis is already associated",ERR,ERROR,*999)
      ELSE
        IF(NUMBER_OF_XI>0.AND.NUMBER_OF_XI<4) THEN
          IF(ANY(XI_DIRECTIONS<1).OR.ANY(XI_DIRECTIONS>3)) CALL FLAG_ERROR("Invalid xi directions specified",ERR,ERROR,*999)
          IF(SIZE(XI_DIRECTIONS,1)/=NUMBER_OF_XI) &
            & CALL FLAG_ERROR("The size of the xi directions array must be the same as the number of xi directions",ERR,ERROR,*999)
          ALLOCATE(NEW_SUB_BASIS,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub-basis",ERR,ERROR,*999)
          NEW_SUB_BASIS%USER_NUMBER=PARENT_BASIS%USER_NUMBER
          NEW_SUB_BASIS%GLOBAL_NUMBER=PARENT_BASIS%GLOBAL_NUMBER
          NEW_SUB_BASIS%FAMILY_NUMBER=PARENT_BASIS%NUMBER_OF_SUB_BASES+1
          NEW_SUB_BASIS%NUMBER_OF_SUB_BASES=0
          NULLIFY(NEW_SUB_BASIS%SUB_BASES)
          NEW_SUB_BASIS%PARENT_BASIS=>PARENT_BASIS
          NEW_SUB_BASIS%NUMBER_OF_XI=NUMBER_OF_XI
          NEW_SUB_BASIS%TYPE=PARENT_BASIS%TYPE
          ALLOCATE(NEW_SUB_BASIS%INTERPOLATION_XI(NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub-basis interpolation xi",ERR,ERROR,*999)
          ALLOCATE(NEW_SUB_BASIS%COLLAPSED_XI(NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub-basis collapsed xi",ERR,ERROR,*999)        
          NUMBER_COLLAPSED=0
          NUMBER_END_COLLAPSED=0
          DO ni=1,NUMBER_OF_XI
            NEW_SUB_BASIS%INTERPOLATION_XI(ni)=PARENT_BASIS%INTERPOLATION_XI(XI_DIRECTIONS(ni))
            NEW_SUB_BASIS%COLLAPSED_XI(ni)=PARENT_BASIS%COLLAPSED_XI(XI_DIRECTIONS(ni))
            IF(NEW_SUB_BASIS%COLLAPSED_XI(ni)==BASIS_XI_COLLAPSED) THEN
              NUMBER_COLLAPSED=NUMBER_COLLAPSED+1
            ELSE IF(NEW_SUB_BASIS%COLLAPSED_XI(ni)==BASIS_COLLAPSED_AT_XI0.OR.NEW_SUB_BASIS%COLLAPSED_XI(ni)== &
              & BASIS_COLLAPSED_AT_XI1) THEN
              NUMBER_END_COLLAPSED=NUMBER_END_COLLAPSED+1
            ENDIF
          ENDDO !ni
          IF(NUMBER_COLLAPSED==0.OR.NUMBER_END_COLLAPSED==0) NEW_SUB_BASIS%COLLAPSED_XI(1:NUMBER_OF_XI)=BASIS_NOT_COLLAPSED
          NULLIFY(NEW_SUB_BASIS%QUADRATURE%BASIS)
          CALL BASIS_QUADRATURE_INITIALISE(NEW_SUB_BASIS,ERR,ERROR,*999)
          NEW_SUB_BASIS%QUADRATURE%TYPE=PARENT_BASIS%QUADRATURE%TYPE
          DO ni=1,NUMBER_OF_XI
            NEW_SUB_BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=PARENT_BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(XI_DIRECTIONS(ni))
          ENDDO !ni
          NEW_SUB_BASIS%BASIS_FINISHED=.TRUE.
          IF(NUMBER_OF_XI>1) THEN
            ALLOCATE(NEW_SUB_BASIS%LINE_BASES(NUMBER_OF_XI),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub-basis line bases",ERR,ERROR,*999)
            IF(NUMBER_OF_XI>2) THEN
              ALLOCATE(NEW_SUB_BASIS%FACE_BASES(NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub-basis face bases",ERR,ERROR,*999)
            ELSE
              ALLOCATE(NEW_SUB_BASIS%FACE_BASES(1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub-basis face bases",ERR,ERROR,*999)
              NEW_SUB_BASIS%FACE_BASES(1)%PTR=>NEW_SUB_BASIS
            ENDIF
          ELSE
            ALLOCATE(NEW_SUB_BASIS%LINE_BASES(1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis line bases",ERR,ERROR,*999)
            NEW_SUB_BASIS%LINE_BASES(1)%PTR=>NEW_SUB_BASIS
            NULLIFY(NEW_SUB_BASIS%FACE_BASES)          
          ENDIF
          !Add the new sub-basis to the list of sub-bases in the parent basis
          ALLOCATE(NEW_SUB_BASES(PARENT_BASIS%NUMBER_OF_SUB_BASES+1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-bases",ERR,ERROR,*999)
          DO basis_idx=1,PARENT_BASIS%NUMBER_OF_SUB_BASES
            NEW_SUB_BASES(basis_idx)%PTR=>PARENT_BASIS%SUB_BASES(basis_idx)%PTR
          ENDDO !basis_idx
          NEW_SUB_BASES(PARENT_BASIS%NUMBER_OF_SUB_BASES+1)%PTR=>NEW_SUB_BASIS
          PARENT_BASIS%NUMBER_OF_SUB_BASES=PARENT_BASIS%NUMBER_OF_SUB_BASES+1
          IF(ASSOCIATED(PARENT_BASIS%SUB_BASES)) DEALLOCATE(PARENT_BASIS%SUB_BASES)
          PARENT_BASIS%SUB_BASES=>NEW_SUB_BASES
          SUB_BASIS=>NEW_SUB_BASIS
        ELSE
          LOCAL_ERROR="Invalid number of xi directions specified ("// &
            & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_XI,"*",ERR,ERROR))//"). You must specify between 1 and 3 xi directions"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Parent basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_SUB_BASIS_CREATE")
    RETURN
999 CALL ERRORS("BASIS_SUB_BASIS_CREATE",ERR,ERROR)
    CALL EXITS("BASIS_SUB_BASIS_CREATE")
    RETURN 1
  END SUBROUTINE BASIS_SUB_BASIS_CREATE
  
  !
  !================================================================================================================================
  !
  
  !>get the type for a basis is identified by a a pointer. \see OPENCMISS::CMISSBasisTypeGet
  SUBROUTINE BASIS_TYPE_GET(BASIS,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the type of the specified basis. \see BASIS_ROUTINES_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BASIS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        TYPE=BASIS%TYPE
      ELSE
        CALL FLAG_ERROR("Basis has not been finished yet",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_TYPE_GET")
    RETURN
999 CALL ERRORS("BASIS_TYPE_GET",ERR,ERROR)
    CALL EXITS("BASIS_TYPE_GET")
    RETURN
  END SUBROUTINE BASIS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a basis is identified by a user number.
  SUBROUTINE BASIS_TYPE_SET_NUMBER(USER_NUMBER,TYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to set.
    INTEGER(INTG), INTENT(IN) :: TYPE !<The type of the basis to set \see BASIS_ROUTINES_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_TYPE_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_TYPE_SET_PTR(BASIS,TYPE,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_TYPE_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_TYPE_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_TYPE_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a basis is identified by a a pointer. \see OPENCMISS::CMISSBasisTypeGet
  SUBROUTINE BASIS_TYPE_SET_PTR(BASIS,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(IN) :: TYPE !<The type of the basis to be set. \see BASIS_ROUTINES_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_TYPE_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          BASIS%TYPE=BASIS_LAGRANGE_HERMITE_TP_TYPE
        CASE(BASIS_SIMPLEX_TYPE)
          !Reset the quadrature
          CALL BASIS_QUADRATURE_FINALISE(BASIS,ERR,ERROR,*999)
          !Change the default parameters for the old basis
          BASIS%TYPE=BASIS_SIMPLEX_TYPE
          BASIS%INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)=BASIS_LINEAR_SIMPLEX_INTERPOLATION
          NULLIFY(BASIS%QUADRATURE%BASIS)
          CALL BASIS_QUADRATURE_INITIALISE(BASIS,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Basis type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is invalid or not implemented"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_TYPE_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_TYPE_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_TYPE_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_TYPE_SET_PTR

  !
  !================================================================================================================================
  !
  
  !>Gets the collapsed xi flags for a basis is identified by a a pointer. \see OPENCMISS::CMISSBasisCollapsedXiGet
  SUBROUTINE BASIS_COLLAPSED_XI_GET(BASIS,COLLAPSED_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). On return, the collapse parameter for each Xi direction. \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_COLLAPSED_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(SIZE(COLLAPSED_XI,1)>=SIZE(BASIS%COLLAPSED_XI)) THEN
          COLLAPSED_XI=BASIS%COLLAPSED_XI
        ELSE
          LOCAL_ERROR="The size of COLLAPSED_XI is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(COLLAPSED_XI,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(BASIS%COLLAPSED_XI,1),"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Basis has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_COLLAPSED_XI_GET")
    RETURN
999 CALL ERRORS("BASIS_COLLAPSED_XI_GET",ERR,ERROR)
    CALL EXITS("BASIS_COLLAPSED_XI_GET")
    RETURN
  END SUBROUTINE BASIS_COLLAPSED_XI_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the collapsed xi flags for a basis is identified by a user number.
  SUBROUTINE BASIS_COLLAPSED_XI_SET_NUMBER(USER_NUMBER,COLLAPSED_XI,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to be set
    INTEGER(INTG), INTENT(IN) :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). The collapse parameter for each Xi direction. \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("BASIS_COLLAPSED_XI_SET_NUMBER",ERR,ERROR,*999)

    CALL BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_COLLAPSED_XI_SET_PTR(BASIS,COLLAPSED_XI,ERR,ERROR,*999)
    
    CALL EXITS("BASIS_COLLAPSED_XI_SET_NUMBER")
    RETURN
999 CALL ERRORS("BASIS_COLLAPSED_XI_SET_NUMBER",ERR,ERROR)
    CALL EXITS("BASIS_COLLAPSED_XI_SET_NUMBER")
    RETURN 1
  END SUBROUTINE BASIS_COLLAPSED_XI_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the collapsed xi flags for a basis is identified by a a pointer. \see OPENCMISS::CMISSBasisCollapsedXiSet
  SUBROUTINE BASIS_COLLAPSED_XI_SET_PTR(BASIS,COLLAPSED_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). The collapse parameter for each Xi direction. \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni1,ni2,ni3,NUMBER_COLLAPSED,COLLAPSED_XI_DIR(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BASIS_COLLAPSED_XI_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FLAG_ERROR("Basis has been finished",ERR,ERROR,*999)
      ELSE
        IF(BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
          IF(BASIS%NUMBER_OF_XI>1) THEN
            IF(SIZE(COLLAPSED_XI,1)==BASIS%NUMBER_OF_XI) THEN
              NUMBER_COLLAPSED=0
              DO ni1=1,BASIS%NUMBER_OF_XI
                SELECT CASE(COLLAPSED_XI(ni1))
                CASE(BASIS_XI_COLLAPSED)
                  NUMBER_COLLAPSED=NUMBER_COLLAPSED+1
                  COLLAPSED_XI_DIR(NUMBER_COLLAPSED)=ni1
                CASE(BASIS_COLLAPSED_AT_XI0,BASIS_COLLAPSED_AT_XI1,BASIS_NOT_COLLAPSED)
                  !Do nothing
                CASE DEFAULT
                  LOCAL_ERROR="Collapsed xi value "//TRIM(NUMBER_TO_VSTRING(COLLAPSED_XI(ni1),"*",ERR,ERROR))// &
                    & " in xi direction "//TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" is invalid"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDDO !ni1
              IF(NUMBER_COLLAPSED>0) THEN
                IF(NUMBER_COLLAPSED<BASIS%NUMBER_OF_XI) THEN
                  IF(BASIS%NUMBER_OF_XI==2) THEN
                    !Two dimensional collapsed basis
                    ni1=COLLAPSED_XI_DIR(1)
                    ni2=OTHER_XI_DIRECTIONS2(ni1)
                    IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI0) THEN
                      IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                        & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                    ELSE IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI1) THEN
                      IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                        & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                    ELSE
                      LOCAL_ERROR="Invalid collapsing of a two dimensional basis. Xi direction "// &
                        & TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" is collapsed so xi direction "// &
                        & TRIM(NUMBER_TO_VSTRING(ni2,"*",ERR,ERROR))//" must be collapsed at an end"
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    !Three dimensional collapsed basis
                    IF(NUMBER_COLLAPSED==1) THEN
                      !One collapse - wedge element
                      ni1=COLLAPSED_XI_DIR(1)
                      ni2=OTHER_XI_DIRECTIONS3(ni1,2,1)
                      ni3=OTHER_XI_DIRECTIONS3(ni1,3,1)
                      IF(COLLAPSED_XI(ni2)==BASIS_NOT_COLLAPSED) THEN
                        IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI0) THEN
                          IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                        ELSE IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI1) THEN
                          IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                        ELSE
                          LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi direction "// &
                            & TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" is collapsed and xi direction "// &
                            & TRIM(NUMBER_TO_VSTRING(ni2,"*",ERR,ERROR))//" is not collapsed so xi direction "// &
                            & TRIM(NUMBER_TO_VSTRING(ni3,"*",ERR,ERROR))//" must be collapsed at an end"
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE IF(COLLAPSED_XI(ni3)==BASIS_NOT_COLLAPSED) THEN
                        IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI0) THEN
                          IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                        ELSE IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI1) THEN
                          IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                        ELSE
                          LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi direction "// &
                            & TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" is collapsed and xi direction "// &
                            & TRIM(NUMBER_TO_VSTRING(ni3,"*",ERR,ERROR))//" is not collapsed so xi direction "// &
                            & TRIM(NUMBER_TO_VSTRING(ni2,"*",ERR,ERROR))//" must be collapsed at an end"
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi direction "// &
                          & TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" is collapsed so one of xi directions "// &
                          & TRIM(NUMBER_TO_VSTRING(ni2,"*",ERR,ERROR))//" or "// &
                          & TRIM(NUMBER_TO_VSTRING(ni3,"*",ERR,ERROR))//" must be collapsed at an end"
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      !Two collapses - pyramid element
                      ni1=COLLAPSED_XI_DIR(1)
                      ni2=COLLAPSED_XI_DIR(2)
                      ni3=OTHER_XI_DIRECTIONS3(ni1,ni2,2)
                      IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI0) THEN
                        IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                          & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                      ELSE IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI1) THEN
                        IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                          & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                      ELSE
                        LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi directions "// &
                          & TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" and "// &
                          & TRIM(NUMBER_TO_VSTRING(ni2,"*",ERR,ERROR))//" are collapsed so xi direction "// &
                          & TRIM(NUMBER_TO_VSTRING(ni3,"*",ERR,ERROR))//" must be collapsed at an end"
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  LOCAL_ERROR="Invalid collapsing of basis. The number of collapsed directions ("// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_COLLAPSED,"*",ERR,ERROR))// &
                    & ") must be less than the number of xi directions ("// &
                    & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//")"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                !No collapses in any xi direction - Reset interpolation_xi if necessary
                DO ni1=1,BASIS%NUMBER_OF_XI
                  IF(BASIS%INTERPOLATION_XI(ni1)==BASIS_QUADRATIC1_HERMITE_INTERPOLATION.OR. &
                    & BASIS%INTERPOLATION_XI(ni1)==BASIS_QUADRATIC2_HERMITE_INTERPOLATION) THEN
                    BASIS%INTERPOLATION_XI(ni1)=BASIS_CUBIC_HERMITE_INTERPOLATION                  
                  ENDIF
                ENDDO
              ENDIF
              BASIS%COLLAPSED_XI=COLLAPSED_XI
            ELSE
              LOCAL_ERROR="The size of the xi collapsed array ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(COLLAPSED_XI,1),"*",ERR,ERROR))//") does not match the number of xi directions ("// &
                & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//") for basis number "// &
                & TRIM(NUMBER_TO_VSTRING(BASIS%USER_NUMBER,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE          
            CALL FLAG_ERROR("Can not collapse a basis with only 1 xi direction",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Can only set collapsed xi directions for a Lagrange Hermite tensor product basis type",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BASIS_COLLAPSED_XI_SET_PTR")
    RETURN
999 CALL ERRORS("BASIS_COLLAPSED_XI_SET_PTR",ERR,ERROR)
    CALL EXITS("BASIS_COLLAPSED_XI_SET_PTR")
    RETURN 1
  END SUBROUTINE BASIS_COLLAPSED_XI_SET_PTR

  !
  !================================================================================================================================
  !

  !>Finds and returns in BASIS a pointer to the basis with the number given in USER_NUMBER. If no basis with that number
  !>exits BASIS is left nullified.
  SUBROUTINE BASIS_USER_NUMBER_FIND(USER_NUMBER,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to find.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On exit, a pointer to the found basis. If no basis with the given user number exists the pointer is NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BASIS_USER_NUMBER_FIND",ERR,ERROR,*999)
    
    CALL BASIS_FAMILY_NUMBER_FIND(USER_NUMBER,0,BASIS,ERR,ERROR,*999)
  
    CALL EXITS("BASIS_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("BASIS_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("BASIS_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE BASIS_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>This routine calculates the weights and abscissae for a Gauss-Legendre quadrature scheme.
  !>\todo Fix this.
  SUBROUTINE GAUSS_LEGENDRE(N,ALPHA,BETA,X,W,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N !<The number of of Gauss points required.
    REAL(DP), INTENT(IN) :: ALPHA !<The lower limit of the integration scheme
    REAL(DP), INTENT(IN) :: BETA !<The upper limit of the integration scheme
    REAL(DP), INTENT(OUT) :: X(N) !<X(i). On exit the i'th Gauss point location
    REAL(DP), INTENT(OUT) :: W(N) !<W(i). On exit the i'th Gauss point weight.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    REAL(DP) :: DIFFERENCE,T1,T2

    INTEGER(INTG) :: GAUSS_START(4) = (/ 0,1,3,6 /)
    REAL(DP) :: XIG(10),WIG(10)
 
!     XIG = (/ 0.500000000000000_DP, &
!       &      0.211324865405187_DP,0.788675134594813_DP, &
!       &      0.112701665379258_DP,0.500000000000000_DP,0.887298334620742_DP, &
!       &      0.06943184420297349_DP,0.330009478207572_DP,0.669990521792428_DP,0.930568155797026_DP /)
!     WIG = (/ 1.000000000000000_DP, &
!       &      0.500000000000000_DP,0.500000000000000_DP, &
!       &      0.277777777777778_DP,0.444444444444444_DP,0.277777777777778_DP, &
!       &      0.173927422568727_DP,0.326072577431273_DP,0.326072577431273_DP,0.173927422568727_DP /)

    XIG = (/ 0.500000000000000_DP, &
      &      (-1.0_DP/sqrt(3.0_DP)+1.0_DP)/2.0_DP,(+1.0_DP/sqrt(3.0_DP)+1.0_DP)/2.0_DP, &
      &      (-SQRT(0.6_DP)+1.0_DP)/2.0_DP, 0.5_DP, (+SQRT(0.6_DP)+1.0_DP)/2.0_DP, &
      &      (-SQRT((3.0_DP+2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP, &      
      &      (-SQRT((3.0_DP-2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP, &
      &      (+SQRT((3.0_DP-2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP, &
      &      (+SQRT((3.0_DP+2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP /)
    WIG = (/ 1.000000000000000_DP, &
      &      0.500000000000000_DP,0.500000000000000_DP, &
      &      2.5_DP/9.0_DP, 4.0_DP/9.0_DP, 2.5_DP/9.0_DP, &
      &      (18.0_DP-SQRT(30.0_DP))/72.0_DP, &
      &      (18.0_DP+SQRT(30.0_DP))/72.0_DP, &
      &      (18.0_DP+SQRT(30.0_DP))/72.0_DP, &
      &      (18.0_DP-SQRT(30.0_DP))/72.0_DP /)
             
    
    CALL ENTERS("GAUSS_LEGENDRE",ERR,ERROR,*999)

    IF(N>1.AND.N<=4) THEN
      DO i=1,N
        X(i)=XIG(GAUSS_START(N)+i)
        W(i)=WIG(GAUSS_START(N)+i)
      ENDDO !i
    ELSE
      CALL FLAG_ERROR("Invalid number of Gauss points. Not implemented",ERR,ERROR,*999)
    ENDIF      
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of gauss points = ",N,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,N,5,5,X,'("Gauss point locations :",5(X,F13.5))','(23X,5(X,F13.5))', &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,N,5,5,W,'("Gauss point weights   :",5(X,F13.5))','(23X,5(X,F13.5))', &
        & ERR,ERROR,*999)
      IF(DIAGNOSTICS2) THEN
        !Check by integrating y=x+1
        T1=0.0_DP !Numerical
        T2=0.0_DP !Analytic
        DO i=1,N
          T1=T1+((X(i)+1.0_DP)*W(i))
        ENDDO !i
        T2=(BETA**2.0_DP/2.0_DP+BETA)-(ALPHA**2.0_DP/2.0_DP-ALPHA)
        DIFFERENCE=ABS(T2-T1)
        CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Numerical Integration Test Difference: ",DIFFERENCE,'(F14.6)', &
          & ERR,ERROR,*999)
      ENDIF
    ENDIF

    CALL EXITS("GAUSS_LEGENDRE")
    RETURN
999 CALL ERRORS("GAUSS_LEGENDRE",ERR,ERROR)
    CALL EXITS("GAUSS_LEGENDRE")
    RETURN 1
  END SUBROUTINE GAUSS_LEGENDRE
  
  !
  !================================================================================================================================
  !

  !>This routine calculates the weights and abscissae for a Gauss quadrature scheme for simplex elements.
  !>
  !>Reference: Liu, Yen and Vinokur, Marcel. "Exact Integrations of Polynomials and Symmetric Quadrature Formulas
  !> over Arbitrary Polyhedral Grids", Journal of Computational Physics, 140:122-147 (1998).
  !>
  SUBROUTINE GAUSS_SIMPLEX(ORDER,NUMBER_OF_VERTICES,N,X,W,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ORDER !<The desired order of the scheme. Currently, the maximum order is 5.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_VERTICES !<The number of vertices. 2, 3 or 4 for lines, triangles or tetrahedra.
    INTEGER(INTG), INTENT(OUT) :: N !<On exit, the number of Gauss points 
    REAL(DP), INTENT(OUT) :: X(:,:) !<X(nj,ng). On exit, the returned positions in area coordinates.
    REAL(DP), INTENT(OUT) :: W(:) !<W(ng). On exit, the returned weights.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ng
    REAL(DP) :: ALPHA_1,ALPHA_2,BETA,LAMBDA,L_C,L1_ALPHA_1,L2_ALPHA_1,L3_ALPHA_1,L4_ALPHA_1,L1_ALPHA_2,L2_ALPHA_2,L3_ALPHA_2, &
      & L4_ALPHA_2,L1_BETA,L2_BETA,L3_BETA,L4_BETA,W_C,W_ALPHA_1,W_ALPHA_2,W_BETA,ACOS_ARG
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GAUSS_SIMPLEX",ERR,ERROR,*999)

    IF(SIZE(X,1)>=NUMBER_OF_VERTICES) THEN
      SELECT CASE(NUMBER_OF_VERTICES)
      CASE(2)
        !Line
        SELECT CASE(ORDER)
        CASE(1)
          N=1
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=1.0_DP
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              W(1)=W_C
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(2)
          N=2
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              ALPHA_1=1.0_DP/SQRT(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP)
              W_ALPHA_1=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=1.0_DP-L1_ALPHA_1
              !Gauss point 1
              X(1,1)=L1_ALPHA_1
              X(2,1)=L2_ALPHA_1
              W(1)=W_ALPHA_1
              !Gauss point 2
              X(1,2)=L2_ALPHA_1
              X(2,2)=L1_ALPHA_1
              W(2)=W_ALPHA_1
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(3)
          N=3
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=-1.0_DP*REAL(NUMBER_OF_VERTICES,DP)*REAL(NUMBER_OF_VERTICES,DP)/ &
                & (REAL(NUMBER_OF_VERTICES,DP)*(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP))
              ALPHA_1=2.0_DP/(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)
              W_ALPHA_1=(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)*(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)/ &
                & (4.0_DP*REAL(NUMBER_OF_VERTICES,DP)*(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP))
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=1.0_DP-L1_ALPHA_1
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              W(1)=W_C
              !Gauss point 2
              X(1,2)=L1_ALPHA_1
              X(2,2)=L2_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L2_ALPHA_1
              X(2,3)=L1_ALPHA_1
              W(3)=W_ALPHA_1
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(4)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(5)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(ORDER,"*",ERR,ERROR))// &
            & " is an invalid Gauss order. You must have an order between 1 and 5"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(3)
        !Triangle
        SELECT CASE(ORDER)
        CASE(1)
          N=1
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=1.0_DP
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              X(3,1)=L_C
              W(1)=W_C
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(2)
          N=3
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              ALPHA_1=-1.0_DP/SQRT(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP)
              W_ALPHA_1=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1
              !Gauss point 1
              X(1,1)=L1_ALPHA_1
              X(2,1)=L2_ALPHA_1
              X(3,1)=L3_ALPHA_1
              W(1)=W_ALPHA_1
              !Gauss point 2
              X(1,2)=L3_ALPHA_1
              X(2,2)=L1_ALPHA_1
              X(3,2)=L2_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L2_ALPHA_1
              X(2,3)=L3_ALPHA_1
              X(3,3)=L1_ALPHA_1
              W(3)=W_ALPHA_1
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(3)
          N=4
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=-1.0_DP*REAL(NUMBER_OF_VERTICES,DP)*REAL(NUMBER_OF_VERTICES,DP)/ &
                & (REAL(NUMBER_OF_VERTICES,DP)*(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP))
              ALPHA_1=2.0_DP/(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)
              W_ALPHA_1=(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)*(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)/ &
                & (4.0_DP*REAL(NUMBER_OF_VERTICES,DP)*(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP))
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              X(3,1)=L_C
              W(1)=W_C
              !Gauss point 2
              X(1,2)=L1_ALPHA_1
              X(2,2)=L2_ALPHA_1
              X(3,2)=L3_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L3_ALPHA_1
              X(2,3)=L1_ALPHA_1
              X(3,3)=L2_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L2_ALPHA_1
              X(2,4)=L3_ALPHA_1
              X(3,4)=L1_ALPHA_1
              W(4)=W_ALPHA_1
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(4)
          N=6
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              ALPHA_1=(-10.0_DP+5.0_DP*SQRT(10.0_DP)+SQRT(950.0_DP-220.0_DP*SQRT(10.0_DP)))/30.0_DP
              ALPHA_2=(-10.0_DP+5.0_DP*SQRT(10.0_DP)-SQRT(950.0_DP-220.0_DP*SQRT(10.0_DP)))/30.0_DP
              W_ALPHA_1=(5.0_DP*ALPHA_2-2.0_DP)/(60.0_DP*ALPHA_1*ALPHA_1*(ALPHA_2-ALPHA_1))
              W_ALPHA_2=(5.0_DP*ALPHA_1-2.0_DP)/(60.0_DP*ALPHA_2*ALPHA_2*(ALPHA_1-ALPHA_2))
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1
              L1_ALPHA_2=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_2=(1.0_DP-ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_2=1.0_DP-L1_ALPHA_2-L2_ALPHA_2
              !Gauss point 1
              X(1,1)=L1_ALPHA_1
              X(2,1)=L2_ALPHA_1
              X(3,1)=L3_ALPHA_1
              W(1)=W_ALPHA_1
              !Gauss point 2
              X(1,2)=L3_ALPHA_1
              X(2,2)=L1_ALPHA_1
              X(3,2)=L2_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L2_ALPHA_1
              X(2,3)=L3_ALPHA_1
              X(3,3)=L1_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L1_ALPHA_2
              X(2,4)=L2_ALPHA_2
              X(3,4)=L3_ALPHA_2
              W(4)=W_ALPHA_2
              !Gauss point 5
              X(1,5)=L3_ALPHA_2
              X(2,5)=L1_ALPHA_2
              X(3,5)=L2_ALPHA_2
              W(5)=W_ALPHA_2
              !Gauss point 6
              X(1,6)=L2_ALPHA_2
              X(2,6)=L3_ALPHA_2
              X(3,6)=L1_ALPHA_2
              W(6)=W_ALPHA_2
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(5)
          N=7
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=9.0_DP/40.0_DP
              ALPHA_1=(1.0_DP+SQRT(15.0_DP))/7.0_DP
              ALPHA_2=(1.0_DP-SQRT(15.0_DP))/7.0_DP
              W_ALPHA_1=(155.0_DP-SQRT(15.0_DP))/1200.0_DP
              W_ALPHA_2=(155.0_DP+SQRT(15.0_DP))/1200.0_DP
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1-L3_ALPHA_1
              L1_ALPHA_2=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_2=(1.0_DP-ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_2=(1.0_DP-ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_2=1.0_DP-L1_ALPHA_2-L2_ALPHA_2-L3_ALPHA_2
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              X(3,1)=L_C
              W(1)=W_C
              !Gauss point 2
              X(1,2)=L1_ALPHA_1
              X(2,2)=L2_ALPHA_1
              X(3,2)=L3_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L3_ALPHA_1
              X(2,3)=L1_ALPHA_1
              X(3,3)=L2_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L2_ALPHA_1
              X(2,4)=L3_ALPHA_1
              X(3,4)=L1_ALPHA_1
              W(4)=W_ALPHA_1
              !Gauss point 5
              X(1,5)=L1_ALPHA_2
              X(2,5)=L2_ALPHA_2
              X(3,5)=L3_ALPHA_2
              W(5)=W_ALPHA_2
              !Gauss point 6
              X(1,6)=L3_ALPHA_2
              X(2,6)=L1_ALPHA_2
              X(3,6)=L2_ALPHA_2
              W(6)=W_ALPHA_2
              !Gauss point 7
              X(1,7)=L2_ALPHA_2
              X(2,7)=L3_ALPHA_2
              X(3,7)=L1_ALPHA_2
              W(7)=W_ALPHA_2
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(ORDER,"*",ERR,ERROR))// &
            & " is an invalid Gauss order. You must have an order between 1 and 5"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(4)
        !Tetrahedra
        SELECT CASE(ORDER)
        CASE(1)
          N=1
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=1.0_DP
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              X(3,1)=L_C
              W(1)=W_C
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(2)
          N=4
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              ALPHA_1=1.0_DP/SQRT(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP)
              W_ALPHA_1=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1-L3_ALPHA_1
              !Gauss point 1
              X(1,1)=L1_ALPHA_1
              X(2,1)=L2_ALPHA_1
              X(3,1)=L3_ALPHA_1
              X(4,1)=L4_ALPHA_1
              W(1)=W_ALPHA_1
              !Gauss point 2
              X(1,2)=L4_ALPHA_1
              X(2,2)=L1_ALPHA_1
              X(3,2)=L2_ALPHA_1
              X(4,2)=L3_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L3_ALPHA_1
              X(2,3)=L4_ALPHA_1
              X(3,3)=L1_ALPHA_1
              X(4,3)=L2_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L2_ALPHA_1
              X(2,4)=L3_ALPHA_1
              X(3,4)=L4_ALPHA_1
              X(4,4)=L1_ALPHA_1
              W(4)=W_ALPHA_1
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(3)
          N=5
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=-1.0_DP*REAL(NUMBER_OF_VERTICES,DP)*REAL(NUMBER_OF_VERTICES,DP)/ &
                & (REAL(NUMBER_OF_VERTICES,DP)*(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP))
              ALPHA_1=2.0_DP/(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)
              W_ALPHA_1=(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)*(REAL(NUMBER_OF_VERTICES,DP)+2.0_DP)/ &
                & (4.0_DP*REAL(NUMBER_OF_VERTICES,DP)*(REAL(NUMBER_OF_VERTICES,DP)+1.0_DP))
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1-L3_ALPHA_1
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              X(3,1)=L_C
              X(4,1)=L_C
              W(1)=W_C
              !Gauss point 2
              X(1,2)=L1_ALPHA_1
              X(2,2)=L2_ALPHA_1
              X(3,2)=L3_ALPHA_1
              X(4,2)=L4_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L4_ALPHA_1
              X(2,3)=L1_ALPHA_1
              X(3,3)=L2_ALPHA_1
              X(4,3)=L3_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L3_ALPHA_1
              X(2,4)=L4_ALPHA_1
              X(3,4)=L1_ALPHA_1
              X(4,4)=L2_ALPHA_1
              W(4)=W_ALPHA_1
              !Gauss point 5
              X(1,5)=L2_ALPHA_1
              X(2,5)=L3_ALPHA_1
              X(3,5)=L4_ALPHA_1
              X(4,5)=L1_ALPHA_1
              W(5)=W_ALPHA_1
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(4)
          N=11
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              L_C=1.0_DP/REAL(NUMBER_OF_VERTICES,DP)
              W_C=-148.0_DP/1875.0_DP
              ALPHA_1=5.0_DP/7.0_DP
              BETA=SQRT(70.0_DP)/28.0_DP
              W_ALPHA_1=343.0_DP/7500.0_DP
              W_BETA=56.0_DP/375.0_DP
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1-L3_ALPHA_1
              L1_BETA=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-2.0_DP)*BETA)/REAL(NUMBER_OF_VERTICES,DP)
              L2_BETA=L1_BETA
              L3_BETA=(1.0_DP-2.0_DP*BETA)/REAL(NUMBER_OF_VERTICES,DP)
              L4_BETA=1.0_DP-L1_BETA-L2_BETA-L3_BETA
              !Gauss point 1
              X(1,1)=L_C
              X(2,1)=L_C
              X(3,1)=L_C
              X(4,1)=L_C
              W(1)=W_C
              !Gauss point 2
              X(1,2)=L1_ALPHA_1
              X(2,2)=L2_ALPHA_1
              X(3,2)=L3_ALPHA_1
              X(4,2)=L4_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L4_ALPHA_1
              X(2,3)=L1_ALPHA_1
              X(3,3)=L2_ALPHA_1
              X(4,3)=L3_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L3_ALPHA_1
              X(2,4)=L4_ALPHA_1
              X(3,4)=L1_ALPHA_1
              X(4,4)=L2_ALPHA_1
              W(4)=W_ALPHA_1
              !Gauss point 5
              X(1,5)=L2_ALPHA_1
              X(2,5)=L3_ALPHA_1
              X(3,5)=L4_ALPHA_1
              X(4,5)=L1_ALPHA_1
              W(5)=W_ALPHA_1
              !Gauss point 6
              X(1,6)=L1_BETA
              X(2,6)=L2_BETA
              X(3,6)=L3_BETA
              X(4,6)=L4_BETA
              W(6)=W_BETA
              !Gauss point 7
              X(1,7)=L1_BETA
              X(2,7)=L3_BETA
              X(3,7)=L2_BETA
              X(4,7)=L4_BETA
              W(7)=W_BETA
              !Gauss point 8
              X(1,8)=L1_BETA
              X(2,8)=L3_BETA
              X(3,8)=L4_BETA
              X(4,8)=L2_BETA
              W(8)=W_BETA
              !Gauss point 9
              X(1,9)=L3_BETA
              X(2,9)=L1_BETA
              X(3,9)=L2_BETA
              X(4,9)=L4_BETA
              W(9)=W_BETA
              !Gauss point 10
              X(1,10)=L3_BETA
              X(2,10)=L1_BETA
              X(3,10)=L4_BETA
              X(4,10)=L2_BETA
              W(10)=W_BETA
              !Gauss point 11
              X(1,11)=L3_BETA
              X(2,11)=L4_BETA
              X(3,11)=L1_BETA
              X(4,11)=L2_BETA
              W(11)=W_BETA
            ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(5)
          N=14
          IF(SIZE(X,2)>=N) THEN
            IF(SIZE(W,1)>=N) THEN
              ACOS_ARG=67.0_DP*SQRT(79.0_DP)/24964.0_DP
              !!todo CHECK THIS!!!
              LAMBDA=4.0_DP/27.0_DP*(4.0_DP*SQRT(79.0_DP)*COS(((ACOS(ACOS_ARG)+TWOPI)/3.0_DP))+71.0_DP)

              ALPHA_1=(SQRT(9.0_DP*LAMBDA*LAMBDA-248.0_DP*LAMBDA+1680.0_DP)+28.0_DP-3.0_DP*LAMBDA)/ &
                & (112.0_DP-10.0_DP*LAMBDA)
              ALPHA_2=(-1.0_DP*SQRT(9.0_DP*LAMBDA*LAMBDA-248.0_DP*LAMBDA+1680.0_DP)+28.0_DP-3.0_DP*LAMBDA)/ &
                & (112.0_DP-10.0_DP*LAMBDA)
              BETA=1.0_DP/SQRT(LAMBDA)
              W_ALPHA_1=((21.0_DP-LAMBDA)*ALPHA_2-7.0_DP)/(420.0_DP*ALPHA_1*ALPHA_1*(ALPHA_2-ALPHA_1))
              W_ALPHA_2=((21.0_DP-LAMBDA)*ALPHA_1-7.0_DP)/(420.0_DP*ALPHA_2*ALPHA_2*(ALPHA_1-ALPHA_2))
              W_BETA=LAMBDA*LAMBDA/840.0_DP
              L1_ALPHA_1=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_1=(1.0_DP-ALPHA_1)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_1=1.0_DP-L1_ALPHA_1-L2_ALPHA_1-L3_ALPHA_1
              L1_ALPHA_2=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-1.0_DP)*ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L2_ALPHA_2=(1.0_DP-ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L3_ALPHA_2=(1.0_DP-ALPHA_2)/REAL(NUMBER_OF_VERTICES,DP)
              L4_ALPHA_2=1.0_DP-L1_ALPHA_2-L2_ALPHA_2-L3_ALPHA_2
              L1_BETA=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-2.0_DP)*BETA)/REAL(NUMBER_OF_VERTICES,DP)
              L2_BETA=(1.0_DP+(REAL(NUMBER_OF_VERTICES,DP)-2.0_DP)*BETA)/REAL(NUMBER_OF_VERTICES,DP)
              L3_BETA=(1.0_DP-2.0_DP*BETA)/REAL(NUMBER_OF_VERTICES,DP)
              L4_BETA=1.0_DP-L1_BETA-L2_BETA-L3_BETA
              !Gauss point 1
              X(1,1)=L1_ALPHA_1
              X(2,1)=L2_ALPHA_1
              X(3,1)=L3_ALPHA_1
              X(4,1)=L4_ALPHA_1
              W(1)=W_ALPHA_1
              !Gauss point 2
              X(1,2)=L4_ALPHA_1
              X(2,2)=L1_ALPHA_1
              X(3,2)=L2_ALPHA_1
              X(4,2)=L3_ALPHA_1
              W(2)=W_ALPHA_1
              !Gauss point 3
              X(1,3)=L3_ALPHA_1
              X(2,3)=L4_ALPHA_1
              X(3,3)=L1_ALPHA_1
              X(4,3)=L2_ALPHA_1
              W(3)=W_ALPHA_1
              !Gauss point 4
              X(1,4)=L2_ALPHA_1
              X(2,4)=L3_ALPHA_1
              X(3,4)=L4_ALPHA_1
              X(4,4)=L1_ALPHA_1
              W(4)=W_ALPHA_1
              !Gauss point 5
              X(1,5)=L1_ALPHA_2
              X(2,5)=L2_ALPHA_2
              X(3,5)=L3_ALPHA_2
              X(4,5)=L4_ALPHA_2
              W(5)=W_ALPHA_2
              !Gauss point 6
              X(1,6)=L4_ALPHA_2
              X(2,6)=L1_ALPHA_2
              X(3,6)=L2_ALPHA_2
              X(4,6)=L3_ALPHA_2
              W(6)=W_ALPHA_2
              !Gauss point 7
              X(1,7)=L3_ALPHA_2
              X(2,7)=L4_ALPHA_2
              X(3,7)=L1_ALPHA_2
              X(4,7)=L2_ALPHA_2
              W(7)=W_ALPHA_2
              !Gauss point 8
              X(1,8)=L2_ALPHA_2
              X(2,8)=L3_ALPHA_2
              X(3,8)=L4_ALPHA_2
              X(4,8)=L1_ALPHA_2
              W(8)=W_ALPHA_2
              !Gauss point 9
              X(1,9)=L1_BETA
              X(2,9)=L2_BETA
              X(3,9)=L3_BETA
              X(4,9)=L4_BETA
              W(9)=W_BETA
              !Gauss point 10
              X(1,10)=L1_BETA
              X(2,10)=L3_BETA
              X(3,10)=L2_BETA
              X(4,10)=L4_BETA
              W(10)=W_BETA
              !Gauss point 11
              X(1,11)=L1_BETA
              X(2,11)=L3_BETA
              X(3,11)=L4_BETA
              X(4,11)=L2_BETA
              W(11)=W_BETA
              !Gauss point 12
              X(1,12)=L3_BETA
              X(2,12)=L1_BETA
              X(3,12)=L2_BETA
              X(4,12)=L4_BETA
              W(12)=W_BETA
              !Gauss point 13
              X(1,13)=L3_BETA
              X(2,13)=L1_BETA
              X(3,13)=L4_BETA
              X(4,13)=L2_BETA
              W(13)=W_BETA
              !Gauss point 14
              X(1,14)=L3_BETA
              X(2,14)=L4_BETA
              X(3,14)=L1_BETA
              X(4,14)=L2_BETA
              W(14)=W_BETA
             ELSE
              LOCAL_ERROR="The first dimension of the W array is "//TRIM(NUMBER_TO_VSTRING(SIZE(W,1),"*",ERR,ERROR))// &
                & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The second dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,2),"*",ERR,ERROR))// &
              & " and it must be >="//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(ORDER,"*",ERR,ERROR))// &
            & " is an invalid Gauss order. You must have an order between 1 and 5"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(NUMBER_OF_VERTICES,"*",ERR,ERROR))// &
          & " is an invalid number of vertices. You must have between 2 and 4 vertices"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      LOCAL_ERROR="The first dimension of the X array is "//TRIM(NUMBER_TO_VSTRING(SIZE(X,1),"*",ERR,ERROR))// &
        & " and it must be >= the number of vertices"
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Simplex Gauss quadrature points:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of vertices = ",NUMBER_OF_VERTICES,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Order = ",ORDER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of gauss points = ",N,ERR,ERROR,*999)
      DO ng=1,N
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Gauss point ",ng,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ORDER,4,4,X(:,ng),'("        Location(nic) :",4(X,F13.5))', &
          & '(23X,4(X,F13.5))',ERR,ERROR,*999)
        CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Weight        : ",W(ng),"F13.5",ERR,ERROR,*999)
      ENDDO !ng
      IF(DIAGNOSTICS2) THEN
!!TODO: \todo add in integral check
      ENDIF
    ENDIF

    CALL EXITS("GAUSS_SIMPLEX")
    RETURN
999 CALL ERRORS("GAUSS_SIMPLEX",ERR,ERROR)
    CALL EXITS("GAUSS_SIMPLEX")
    RETURN 1
  END SUBROUTINE GAUSS_SIMPLEX
  
  !
  !================================================================================================================================
  !

  !>Evaluates a 1D cubic Hermite basis function.
  FUNCTION HERMITE_CUBIC_EVALUATE(NODE_INDEX,NODE_DERIVATIVE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,ERR,ERROR)
  
   !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: NODE_DERIVATIVE_INDEX !<The local derivative number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: HERMITE_CUBIC_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    CALL ENTERS("HERMITE_CUBIC_EVALUATE",ERR,ERROR,*999)

    HERMITE_CUBIC_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=(2.0_DP*XI-3.0_DP)*XI*XI+1.0_DP ! 2xi^3-3xi^2+1
        CASE(2)
          HERMITE_CUBIC_EVALUATE=((XI-2.0_DP)*XI+1.0_DP)*XI ! xi^3-2xi^2+xi
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
        END SELECT
      CASE(2)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=XI*XI*(3.0_DP-2.0_DP*XI) ! -2xi^3+3xi^2
        CASE(2)
          HERMITE_CUBIC_EVALUATE=XI*XI*(XI-1.0_DP) ! xi^3-xi^2
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI*(XI-1.0_DP) ! 6xi^2-6xi
        CASE(2)
          HERMITE_CUBIC_EVALUATE=(3.0_DP*XI-4.0_DP)*XI+1.0_DP ! 3xi^2-4xi+1
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
        END SELECT
      CASE(2)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI*(1.0_DP-XI) ! -6xi^2+6xi
        CASE(2)
          HERMITE_CUBIC_EVALUATE=XI*(3.0_DP*XI-2.0_DP) ! 3xi^2-2xi
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=12.0_DP*XI-6.0_DP ! 12xi-6
        CASE(2)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI-4.0_DP ! 6xi-4
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
        END SELECT
      CASE(2)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=6.0_DP-12.0_DP*XI ! -12xi+6
        CASE(2)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI-2.0_DP ! 6xi-2
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("HERMITE_CUBIC_EVALUATE")
    RETURN
999 CALL ERRORS("HERMITE_CUBIC_EVALUATE",ERR,ERROR)
    CALL EXITS("HERMITE_CUBIC_EVALUATE")
    RETURN 
  END FUNCTION HERMITE_CUBIC_EVALUATE
 
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: HERMITE_QUADRATIC_EVALUATE
  !###  Description:
  !###    Evaluates a 1D quadratic Hermite basis function at position XI,and with the give NODE_INDEX, NODE_DERIVATIVE_INDEX and
  !###    PARTIAL_DERIVATIVE_INDEX. SPECIAL_NODE_INDEX is the node with no derivative term.
  !###  Child-functions: HERMITE_QUADRATIC_EVALUATE_DP

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D quadratic Hermite basis function
  FUNCTION HERMITE_QUADRATIC_EVALUATE(NODE_INDEX,NODE_DERIVATIVE_INDEX,PARTIAL_DERIVATIVE_INDEX,SPECIAL_NODE_INDEX,XI,ERR,ERROR)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: NODE_DERIVATIVE_INDEX !<The local derivative number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index to evaluate.
    INTEGER(INTG), INTENT(IN) :: SPECIAL_NODE_INDEX !<The local node number with no derivative term.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: HERMITE_QUADRATIC_EVALUATE
    !Local variables
    
    CALL ENTERS("HERMITE_QUADRATIC_EVALUATE",ERR,ERROR,*999)
    
    HERMITE_QUADRATIC_EVALUATE=0.0_DP
    SELECT CASE(SPECIAL_NODE_INDEX)
    CASE(1)
      SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
      CASE(NO_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=(XI-2.0_DP)*XI+1.0_DP ! xi^2-2xi+1
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=(2.0_DP-XI)*XI ! -xi^2+2xi
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=(XI-1.0_DP)*XI ! xi^2-xi
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
        END SELECT
      CASE(FIRST_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP*XI-2.0_DP  ! 2xi-2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP*XI+2.0_DP ! -2xi+2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP*XI-1.0_DP ! 2xi-1
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
        END SELECT
      CASE(SECOND_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP ! 2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP ! -2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP ! 2
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
      END SELECT
    CASE(2)
      SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
      CASE(NO_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=1.0_DP-XI*XI ! -xi^2+1
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=XI*(1.0_DP-XI) ! -xi^2+xi
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=XI*XI ! xi^2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
        END SELECT
      CASE(FIRST_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP*XI ! -2xi
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=1.0_DP-2.0_DP*XI ! -2xi+1
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP*XI ! 2xi
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
        END SELECT
      CASE(SECOND_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP ! -2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP ! -2
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP ! 2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid node derivative index",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid special node index",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("HERMITE_QUADRATIC_EVALUATE")
    RETURN
999 CALL ERRORS("HERMITE_QUADRATIC_EVALUATE",ERR,ERROR)
    CALL EXITS("HERMITE_QUADRATIC_EVALUATE")
    RETURN 
  END FUNCTION HERMITE_QUADRATIC_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D cubic Lagrange basis function.
  FUNCTION LAGRANGE_CUBIC_EVALUATE(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node of the basis to evaluate. Must be between 1 and 4.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: LAGRANGE_CUBIC_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    CALL ENTERS("LAGRANGE_CUBIC_EVALUATE",ERR,ERROR,*999)
    
    LAGRANGE_CUBIC_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_CUBIC_EVALUATE=0.5_DP*(3.0_DP*XI-1.0_DP)*(3.0_DP*XI-2.0_DP)*(1.0_DP-XI) !
      CASE(2)
        LAGRANGE_CUBIC_EVALUATE=4.5_DP*XI*(3.0_DP*XI-2.0_DP)*(XI-1.0_DP) !
      CASE(3)
        LAGRANGE_CUBIC_EVALUATE=4.5_DP*XI*(3.0_DP*XI-1.0_DP)*(1.0_DP-XI) !
      CASE(4)
        LAGRANGE_CUBIC_EVALUATE=0.5_DP*XI*(3.0_DP*XI-1.0_DP)*(3.0_DP*XI-2.0_DP) !
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_CUBIC_EVALUATE=-13.5_DP*XI*XI+18.0_DP*XI-5.5_DP ! -13.5xi^2+18xi-5.5
      CASE(2)
        LAGRANGE_CUBIC_EVALUATE= 40.5_DP*XI*XI-45.0_DP*XI+9.0_DP ! 40.5xi^2-45xi+9
      CASE(3)
        LAGRANGE_CUBIC_EVALUATE=-40.5_DP*XI*XI+36.0_DP*XI-4.5_DP ! -40.5xi^2+36xi-4.5
      CASE(4)
        LAGRANGE_CUBIC_EVALUATE= 13.5_DP*XI*XI- 9.0_DP*XI+1.0_DP ! 13.5xi^2-9xi+1
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(2.0_DP-3.0_DP*XI) ! 18-27xi
      CASE(2)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(9.0_DP*XI-5.0_DP) ! 81xi-45
      CASE(3)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(4.0_DP-9.0_DP*XI) ! 36-81xi
      CASE(4)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(3.0_DP*XI-1.0_DP) ! 27xi-9
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("LAGRANGE_CUBIC_EVALUATE")
    RETURN
999 CALL ERRORS("LAGRANGE_CUBIC_EVALUATE",ERR,ERROR)
    CALL EXITS("LAGRANGE_CUBIC_EVALUATE")
    RETURN 
  END FUNCTION LAGRANGE_CUBIC_EVALUATE

  !
  !================================================================================================================================
  !
  
  !> Evaluates a 1D linear Lagrange basis function.
  FUNCTION LAGRANGE_LINEAR_EVALUATE(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,ERR,ERROR)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node of the basis to evaluate. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: LAGRANGE_LINEAR_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    CALL ENTERS("LAGRANGE_LINEAR_EVALUATE",ERR,ERROR,*999)

    LAGRANGE_LINEAR_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_LINEAR_EVALUATE=1.0_DP-XI ! 1-xi
      CASE(2)
        LAGRANGE_LINEAR_EVALUATE=XI !xi
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_LINEAR_EVALUATE=-1.0_DP ! -1
      CASE(2)
        LAGRANGE_LINEAR_EVALUATE=1.0_DP ! 1
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_LINEAR_EVALUATE=0.0_DP ! 0
      CASE(2)
        LAGRANGE_LINEAR_EVALUATE=0.0_DP ! 0
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
    END SELECT
    
    CALL EXITS("LAGRANGE_LINEAR_EVALUATE")
    RETURN
999 CALL ERRORS("LAGRANGE_LINEAR_EVALUATE",ERR,ERROR)
    CALL EXITS("LAGRANGE_LINEAR_EVALUATE")
    RETURN 
  END FUNCTION LAGRANGE_LINEAR_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D quadratic Lagrange basis function.
  FUNCTION LAGRANGE_QUADRATIC_EVALUATE(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,ERR,ERROR)
     
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node of the basis to evaluate. Must be between 1 and 3.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: LAGRANGE_QUADRATIC_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    CALL ENTERS("LAGRANGE_QUADRATIC_EVALUATE",ERR,ERROR,*999)

    LAGRANGE_QUADRATIC_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_QUADRATIC_EVALUATE=1.0_DP-3.0_DP*XI+2.0_DP*XI*XI ! 1-3xi+2xi^2
      CASE(2)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP*XI*(1.0_DP-XI) ! 4xi-4xi^2
      CASE(3)
        LAGRANGE_QUADRATIC_EVALUATE=XI*(XI+XI-1.0_DP) ! 2xi^2-xi
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP*XI-3.0_DP ! 4xi-3
      CASE(2)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP-8.0_DP*XI ! 4-8xi
      CASE(3)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP*XI-1.0_DP ! 4xi-1
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP ! 4
      CASE(2)
        LAGRANGE_QUADRATIC_EVALUATE=-8.0_DP ! -8
      CASE(3)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP ! 4
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("LAGRANGE_QUADRATIC_EVALUATE")
    RETURN
999 CALL ERRORS("LAGRANGE_QUADRATIC_EVALUATE",ERR,ERROR)
    CALL EXITS("LAGRANGE_QUADRATIC_EVALUATE")
    RETURN 
  END FUNCTION LAGRANGE_QUADRATIC_EVALUATE

  !
  !================================================================================================================================
  !
  
  !>Evaluates a cubic simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments
  FUNCTION SIMPLEX_CUBIC_EVALUATE_DP(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XL,ERR,ERROR)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: XL !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: SIMPLEX_CUBIC_EVALUATE_DP
    !Local variables
    
    CALL ENTERS("SIMPLEX_CUBIC_EVALUATE_DP",ERR,ERROR,*999)
    
    SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP
        
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=1.0 !1
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP*XL !3L
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP/2.0_DP*XL*(3.0_DP*XL-1.0_DP) !3/2.L(3L-1)
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=0.5_DP*XL*(3.0_DP*XL-1.0_DP)*(3.0_DP*XL-2.0_DP) !1/2.L(3L-1)(3L-2)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP !3
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP/2.0_DP*(6.0_DP*XL-1) !3/2.(6L-1)
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=13.5_DP*XL*XL-9.0_DP*XL+1.0_DP !27/2.L^2-9L+1
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=9.0_DP !9
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=2.0_DP*XL-9.0_DP
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=2.0_DP !2
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index.",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("SIMPLEX_CUBIC_EVALUATE_DP")
    RETURN
999 CALL ERRORS("SIMPLEX_CUBIC_EVALUATE_DP",ERR,ERROR)
    CALL EXITS("SIMPLEX_CUBIC_EVALUATE_DP")
    RETURN 
  END FUNCTION SIMPLEX_CUBIC_EVALUATE_DP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a linear simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments
  FUNCTION SIMPLEX_LINEAR_EVALUATE_DP(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XL,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: XL !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: SIMPLEX_LINEAR_EVALUATE_DP
    !Local variables
    
    CALL ENTERS("SIMPLEX_LINEAR_EVALUATE_DP",ERR,ERROR,*999)

    SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=1.0 !1
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=XL  !L
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP  !0
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=1.0_DP !1
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index",ERR,ERROR,*999)
    END SELECT
    
    CALL EXITS("SIMPLEX_LINEAR_EVALUATE_DP")
    RETURN
999 CALL ERRORS("SIMPLEX_LINEAR_EVALUATE_DP",ERR,ERROR)
    CALL EXITS("SIMPLEX_LINEAR_EVALUATE_DP")
    RETURN 
  END FUNCTION SIMPLEX_LINEAR_EVALUATE_DP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a quadratic simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments.
   FUNCTION SIMPLEX_QUADRATIC_EVALUATE_DP(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XL,ERR,ERROR)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: XL !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: SIMPLEX_QUADRATIC_EVALUATE_DP
    !Local variables
    
    CALL ENTERS("SIMPLEX_QUADRATIC_EVALUATE_DP",ERR,ERROR,*999)

    SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=1.0_DP !1
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=2.0_DP*XL !2L
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=XL*(2.0_DP*XL-1.0_DP) !L(2L-1)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=2.0_DP !4
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=4.0_DP*XL-1.0_DP !4L-1
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index",ERR,ERROR,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=4.0_DP !4
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid node index.",ERR,ERROR,*999)
      END SELECT
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid partial derivative index.",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("SIMPLEX_QUADRATIC_EVALUATE_DP")
    RETURN
999 CALL ERRORS("SIMPLEX_QUADRATIC_EVALUATE_DP",ERR,ERROR)
    CALL EXITS("SIMPLEX_QUADRATIC_EVALUATE_DP")
    RETURN 
  END FUNCTION SIMPLEX_QUADRATIC_EVALUATE_DP

  !
  !================================================================================================================================
  !

END MODULE BASIS_ROUTINES

