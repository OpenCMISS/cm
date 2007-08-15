!> \file
!> $Id: coordinate_routines.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module contains all coordinate transformation and support routines.
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

!> This module contains all coordinate transformation and support routines.
MODULE COORDINATE_ROUTINES

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE MATHS
  USE STRINGS
  USE TYPES
  
  IMPLICIT NONE

!!TODO: Should all of the get/set/create/destroy routines be accessed via pointers???

  PRIVATE

  !Module parameters

  !> \addtogroup COORINDATE_ROUTINES_CoordinateSystemTypes COORDINATE_ROUTINES::CoordinateSystemTypes
  !> \see COORDINATE_ROUTINES
  !> Coordinate system type parameters
  !>@{ 
  INTEGER(INTG), PARAMETER :: COORDINATE_RECTANGULAR_CARTESIAN_TYPE=1 !<Rectangular Cartesian coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_CYCLINDRICAL_POLAR_TYPE=2 !<Cylindrical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_SPHERICAL_POLAR_TYPE=3 !<Spherical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_PROLATE_SPHEROIDAL_TYPE=4 !<Prolate spheroidal coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_OBLATE_SPHEROIDAL_TYPE=5 !<Oblate spheroidal coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  !>@}

  !> \addtogroup COORDINATE_ROUTINES_RadialInterpolations COORDINATE_ROUTINES::RadialInterpolations
  !> \see COORDINATE_ROUTINES
  !> \brief The type of radial interpolation for polar coordinate systems
  !>@{
  INTEGER(INTG), PARAMETER :: COORDINATE_NO_RADIAL_INTERPOLATION_TYPE=0 !<No radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_INTERPOLATION_TYPE=1 !<r radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE=2 !<r^2 radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE=3 !<r^3 radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  !>@}
  
  !> \addtogroup COORDINATE_ROUTINES_JacobianType COORDINATE_ROUTINES::JacobianType
  !> \see COORDINATE_ROUTINES
  !> \brief The type of Jacobian to return when coordinate metrics are calculated.
  !>@{
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_LINE_TYPE=1 !<Line type Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_AREA_TYPE=2 !<Area type Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_VOLUME_TYPE=3 !<Volume type Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  !>@}
  
  !Module types

  TYPE COORDINATE_SYSTEM_PTR_TYPE
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: PTR
  END TYPE COORDINATE_SYSTEM_PTR_TYPE
  
  TYPE COORDINATE_SYSTEMS_TYPE
    INTEGER(INTG) :: NUMBER_OF_COORDINATE_SYSTEMS
    TYPE(COORDINATE_SYSTEM_PTR_TYPE), POINTER :: COORDINATE_SYSTEMS(:)
  END TYPE COORDINATE_SYSTEMS_TYPE
  
  !Module variables

  CHARACTER(LEN=21) :: COORDINATE_SYSTEM_TYPE_STRING(5) = &
    & (/ "Rectangular Cartesian",&
    &    "Cylindrical Polar    ", &
    &    "Spherical Polar      ", &
    &    "Prolate Spheroidal   ", &
    &    "Oblate Spheroidal    " /)

  TYPE(COORDINATE_SYSTEMS_TYPE) :: COORDINATE_SYSTEMS
  
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: GLOBAL_COORDINATE_SYSTEM
  
  !Interfaces

  INTERFACE COORDINATE_CONVERT_FROM_RC
    MODULE PROCEDURE COORDINATE_CONVERT_FROM_RC_DP
    MODULE PROCEDURE COORDINATE_CONVERT_FROM_RC_SP
  END INTERFACE !COORDINATE_CONVERT_FROM_RC

  INTERFACE COORDINATE_CONVERT_TO_RC
    MODULE PROCEDURE COORDINATE_CONVERT_TO_RC_DP
    MODULE PROCEDURE COORDINATE_CONVERT_TO_RC_SP
  END INTERFACE !COORDINATE_CONVERT_TO_RC

  INTERFACE COORDINATE_DELTA_CALCULATE
    MODULE PROCEDURE COORDINATE_DELTA_CALCULATE_DP
    !MODULE PROCEDURE COORDINATE_DELTA_CALCULATE_SP
  END INTERFACE !COORDINATE_DELTA_CALCULATE

  INTERFACE COORDINATE_SYSTEM_DIMENSION_SET
    MODULE PROCEDURE COORDINATE_SYSTEM_DIMENSION_SET_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_DIMENSION_SET_PTR
  END INTERFACE !COORDINATE_SYSTEM_DIMENSION_SET

  INTERFACE COORDINATE_SYSTEM_FOCUS_SET
    MODULE PROCEDURE COORDINATE_SYSTEM_FOCUS_SET_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_FOCUS_SET_PTR
  END INTERFACE !COORDINATE_SYSTEM_FOCUS_SET

  INTERFACE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET
    MODULE PROCEDURE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR
  END INTERFACE !COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET

  INTERFACE COORDINATE_SYSTEM_TYPE_SET
    MODULE PROCEDURE COORDINATE_SYSTEM_TYPE_SET_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_TYPE_SET_PTR
  END INTERFACE !COORDINATE_SYSTEM_TYPE_SET

  INTERFACE COORDINATE_SYSTEM_ORIGIN_SET
    MODULE PROCEDURE COORDINATE_SYSTEM_ORIGIN_SET_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_ORIGIN_SET_PTR
  END INTERFACE !COORDINATE_SYSTEM_ORIGIN_SET

  INTERFACE COORDINATE_SYSTEM_ORIENTATION_SET
    MODULE PROCEDURE COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_ORIENTATION_SET_PTR
  END INTERFACE !COORDINATE_SYSTEM_ORIENTATION_SET

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE DXZ
    MODULE PROCEDURE DXZ_DP
    !MODULE PROCEDURE DXZ_SP
  END INTERFACE !Dxz

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE D2ZX
    MODULE PROCEDURE D2ZX_DP
    !MODULE PROCEDURE D2ZX_SP
  END INTERFACE !D2ZX

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE DZX
    MODULE PROCEDURE DZX_DP
    !MODULE PROCEDURE DZX_SP
  END INTERFACE !DZX

  INTERFACE COORDINATE_DERIVATIVE_CONVERT_TO_RC
    MODULE PROCEDURE COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP
    MODULE PROCEDURE COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP
  END INTERFACE !COORDINATE_DERIVATIVE_CONVERT_TO_RC

  INTERFACE COORDINATE_SYSTEM_DESTROY
    MODULE PROCEDURE COORDINATE_SYSTEM_DESTROY_NUMBER
    MODULE PROCEDURE COORDINATE_SYSTEM_DESTROY_PTR
  END INTERFACE !COORDINATE_SYSTEM_DESTROY

  PUBLIC COORDINATE_RECTANGULAR_CARTESIAN_TYPE,COORDINATE_CYCLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE, &
    & COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE

  PUBLIC COORDINATE_NO_RADIAL_INTERPOLATION_TYPE,COORDINATE_RADIAL_INTERPOLATION_TYPE, &
    & COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE,COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE
  
  PUBLIC COORDINATE_JACOBIAN_LINE_TYPE,COORDINATE_JACOBIAN_AREA_TYPE,COORDINATE_JACOBIAN_VOLUME_TYPE
  
  PUBLIC COORDINATE_SYSTEM_TYPE_STRING

  PUBLIC COORDINATE_CONVERT_FROM_RC,COORDINATE_CONVERT_TO_RC,COORDINATE_DELTA_CALCULATE,COORDINATE_DERIVATIVE_NORM, &
    & COORDINATE_INTERPOLATION_ADJUST,COORDINATE_INTERPOLATION_PARAMETERS_ADJUST,COORDINATE_METRICS_CALCULATE, &
    & COORDINATE_SYSTEM_DIMENSION_GET,COORDINATE_SYSTEM_FOCUS_GET,COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET, &
    & COORDINATE_SYSTEM_TYPE_GET,COORDINATE_SYSTEM_DIMENSION_SET,COORDINATE_SYSTEM_FOCUS_SET, &
    & COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET,COORDINATE_SYSTEM_TYPE_SET,COORDINATE_SYSTEM_ORIGIN_SET, &
    & COORDINATE_SYSTEM_ORIENTATION_SET,COORDINATE_SYSTEM_CREATE_START,COORDINATE_SYSTEM_CREATE_FINISH, &
    & COORDINATE_DERIVATIVE_CONVERT_TO_RC,COORDINATE_SYSTEM_DESTROY,COORDINATE_SYSTEM_USER_NUMBER_FIND, &
    & COORDINATE_SYSTEMS_INITIALISE,COORDINATE_SYSTEMS_FINALISE
  
CONTAINS

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: COORDINATE_CONVERT_FROM_RC
  !###  Description:
  !###    COORDINATE_CONVERT_FROM_RC performs a coordinate transformation from a rectangular cartesian coordinate at the point with
  !###    coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by COORDINATE_SYSTEM.
  !###  Child-functions: COORDINATE_CONVERT_FROM_RC_DP,COORDINATE_CONVERT_FROM_RC_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION COORDINATE_CONVERT_FROM_RC_DP(COORDINATE_SYSTEM,Z,ERR,ERROR)
  
    !#### Function: COORDINATE_CONVERT_FROM_RC_DP
    !###  Type: REAL(DP)(SIZE(Z,1))
    !###  Description:
    !###    COORDINATE_CONVERT_FROM_RC_DP performs a coordinate transformation from a rectangular cartesian coordinate at the
    !###    point with coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by
    !###    COORDINATE_SYSTEM for double precision coordinates.
    !###  Parent-function: COORDINATE_CONVERT_FROM_RC
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    REAL(DP), INTENT(IN) :: Z(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: COORDINATE_CONVERT_FROM_RC_DP(SIZE(Z,1))
    !Local variables
    REAL(DP) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_FROM_RC_DP",ERR,ERROR,*999)

    COORDINATE_CONVERT_FROM_RC_DP=0.0_DP

    IF(SIZE(Z,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of Z is less than the number of dimensions",ERR,ERROR,*999)
    
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_FROM_RC_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=Z(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(1),Z(2))
      CASE(3)
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(1),Z(2))
        COORDINATE_CONVERT_FROM_RC_DP(3)=Z(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
      IF(COORDINATE_CONVERT_FROM_RC_DP(2)<0.0_DP) &
        & COORDINATE_CONVERT_FROM_RC_DP(2)=COORDINATE_CONVERT_FROM_RC_DP(2)+2.0_DP*PI !reference coordinate 0->2*pi
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2+Z(3)**2)
        IF(Z(1)/=0.0_DP.OR.Z(2)/=0.0_DP) THEN
          COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(2),Z(1))
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(2)=0.0_DP
        ENDIF
        A1=SQRT(Z(1)**2+Z(2)**2)
        IF(Z(3)/=0.0_DP.OR.A1/=0.0_DP) THEN
          COORDINATE_CONVERT_FROM_RC_DP(3)=ATAN2(Z(3),A1)
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(3)=0.0_DP
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      FOCUS=COORDINATE_SYSTEM%FOCUS
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        A1=Z(1)**2+Z(2)**2+Z(3)**2-FOCUS**2
        A2=SQRT(A1**2+4.0_DP*(FOCUS**2)*(Z(2)**2+Z(3)**2))
        A3=2.0_DP*FOCUS**2
        A4=MAX((A2+A1)/A3,0.0_DP)
        A5=MAX((A2-A1)/A3,0.0_DP)
        A6=SQRT(A4)
        A7=MIN(SQRT(A5),1.0_DP)
        IF(ABS(A7)<=1.0_DP) THEN
          A8=ASIN(A7)
        ELSE
          A8=0.0_DP
          CALL FLAG_WARNING("Put A8=0 since ABS(A8)>1",ERR,ERROR,*999)
        ENDIF
        IF((Z(3)==0.0_DP).OR.(A6==0.0_DP).OR.(A7==0.0_DP)) THEN
          A9=0.0_DP
        ELSE
          IF(ABS(A6*A7)>0.0_DP) THEN
            A9=Z(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0_DP
            CALL FLAG_WARNING("Put A9=0 since A6*A7=0",ERR,ERROR,*999)
          ENDIF
          IF(A9>=1.0_DP) THEN
            A9=PI/2.0_DP
          ELSE IF(A9<=-1.0_DP) THEN
            A9=-PI/2.0_DP
          ELSE
            A9=ASIN(A9)
          ENDIF
        ENDIF
        COORDINATE_CONVERT_FROM_RC_DP(1)=LOG(A6+SQRT(A4+1.0_DP))
        IF(Z(1)>=0.0_DP) THEN
          COORDINATE_CONVERT_FROM_RC_DP(2)=A8
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(2)=PI-A8
        ENDIF
        IF(Z(2)>=0.0_DP) THEN
          COORDINATE_CONVERT_FROM_RC_DP(3)=MOD(A9+2.0_DP*PI,2.0_DP*PI)
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(3)=PI-A9
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("COORDINATE_CONVERT_FROM_RC_DP")
    RETURN
999 CALL ERRORS("COORDINATE_CONVERT_FROM_RC_DP",ERR,ERROR)
    CALL EXITS("COORDINATE_CONVERT_FROM_RC_DP")
    RETURN 
  END FUNCTION COORDINATE_CONVERT_FROM_RC_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION COORDINATE_CONVERT_FROM_RC_SP(COORDINATE_SYSTEM,Z,ERR,ERROR)
  
    !#### Function: COORDINATE_CONVERT_FROM_RC_SP
    !###  Type: REAL(SP)(SIZE(Z,1))
    !###  Description:
    !###    COORDINATE_CONVERT_FROM_RC_SP performs a coordinate transformation from a rectangular cartesian coordinate at the
    !###    point with coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by
    !###    COORDINATE_SYSTEM for single precision coordinates.
    !###  Parent-function: COORDINATE_CONVERT_FROM_RC
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    REAL(SP), INTENT(IN) :: Z(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(SP) :: COORDINATE_CONVERT_FROM_RC_SP(SIZE(Z,1))
    !Local variables
    REAL(SP) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_FROM_RC_SP",ERR,ERROR,*999)

    COORDINATE_CONVERT_FROM_RC_SP=0.0_SP
    
    IF(SIZE(Z,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of Z is less than the number of dimensions",ERR,ERROR,*999)
    
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_FROM_RC_SP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=Z(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(1),Z(2))
      CASE(3)
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(1),Z(2))
        COORDINATE_CONVERT_FROM_RC_SP(3)=Z(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
      IF(COORDINATE_CONVERT_FROM_RC_SP(2)<0.0_SP)  &
        & COORDINATE_CONVERT_FROM_RC_SP(2)=COORDINATE_CONVERT_FROM_RC_SP(2)+2.0_SP*REAL(PI,SP) !reference coordinate 0->2*pi
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2+Z(3)**2)
        IF(Z(1)/=0.0_SP.OR.Z(2)/=0.0_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(2),Z(1))
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(2)=0.0_SP
        ENDIF
        A1=SQRT(Z(1)**2+Z(2)**2)
        IF(Z(3)/=0.0_SP.OR.A1/=0.0_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(3)=ATAN2(Z(3),A1)
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(3)=0.0_SP
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        A1=Z(1)**2+Z(2)**2+Z(3)**2-FOCUS**2
        A2=SQRT(A1**2+4.0_SP*(FOCUS**2)*(Z(2)**2+Z(3)**2))
        A3=2.0_SP*FOCUS**2
        A4=MAX((A2+A1)/A3,0.0_SP)
        A5=MAX((A2-A1)/A3,0.0_SP)
        A6=SQRT(A4)
        A7=MIN(SQRT(A5),1.0_SP)
        IF(ABS(A7)<=1.0_SP) THEN
          A8=ASIN(A7)
        ELSE
          A8=0.0_SP
          CALL FLAG_WARNING("Put A8=0 since ABS(A8)>1",ERR,ERROR,*999)
        ENDIF
        IF((Z(3)==0.0_SP).OR.(A6==0.0_SP).OR.(A7==0.0_SP)) THEN
          A9=0.0_SP
        ELSE
          IF(ABS(A6*A7)>0.0_SP) THEN
            A9=Z(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0_SP
            CALL FLAG_WARNING("Put A9=0 since A6*A7=0",ERR,ERROR,*999)
          ENDIF
          IF(A9>=1.0_SP) THEN
            A9=REAL(PI,SP)/2.0_SP
          ELSE IF(A9<=-1.0_SP) THEN
            A9=-REAL(PI,SP)/2.0_SP
          ELSE
            A9=ASIN(A9)
          ENDIF
        ENDIF
        COORDINATE_CONVERT_FROM_RC_SP(1)=LOG(A6+SQRT(A4+1.0_SP))
        IF(Z(1)>=0.0_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(2)=A8
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(2)=REAL(PI,SP)-A8
        ENDIF
        IF(Z(2)>=0.0_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(3)=MOD(A9+2.0_SP*REAL(PI,SP),2.0_SP*&
            & REAL(PI,SP))
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(3)=REAL(PI,SP)-A9
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("COORDINATE_CONVERT_FROM_RC_SP")
    RETURN
999 CALL ERRORS("COORDINATE_CONVERT_FROM_RC_SP",ERR,ERROR)
    CALL EXITS("COORDINATE_CONVERT_FROM_RC_SP")
    RETURN 
  END FUNCTION COORDINATE_CONVERT_FROM_RC_SP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: COORDINATE_CONVERT_TO_RC
  !###  Description:
  !###    COORDINATE_CONVERT_TO_RC performs a coordinate transformation from a coordinate system identified by COORDINATE_SYSTEM
  !###    at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates.
  !###  Child-functions: COORDINATE_CONVERT_TO_RC_DP,COORDINATE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION COORDINATE_CONVERT_TO_RC_DP(COORDINATE_SYSTEM,X,ERR,ERROR)
  
    !#### Function: COORDINATE_CONVERT_TO_RC_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    COORDINATE_CONVERT_TO_RC_DP performs a coordinate transformation from a coordinate system identified by
    !###    COORDINATE_SYSTEM at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates for
    !###    double precision coordinates.
    !###  Parent-function: COORDINATE_CONVERT_TO_RC
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: COORDINATE_CONVERT_TO_RC_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_TO_RC_DP",ERR,ERROR,*999)
    
    COORDINATE_CONVERT_TO_RC_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions",ERR,ERROR,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_TO_RC_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=X(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))
      CASE(3)
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_DP(3)=X(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN  
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(3)=X(1)*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        COORDINATE_CONVERT_TO_RC_DP(1)=FOCUS*COSH(X(1))*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        COORDINATE_CONVERT_TO_RC_DP(1)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("COORDINATE_CONVERT_TO_RC_DP")
    RETURN
999 CALL ERRORS("COORDINATE_CONVERT_TO_RC_DP",ERR,ERROR)
    CALL EXITS("COORDINATE_CONVERT_TO_RC_DP")
    RETURN 
  END FUNCTION COORDINATE_CONVERT_TO_RC_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION COORDINATE_CONVERT_TO_RC_SP(COORDINATE_SYSTEM,X,ERR,ERROR)
  
    !#### Function: COORDINATE_CONVERT_TO_RC_SP
    !###  Type: REAL(SP)(SIZE(X,1))
    !###  Description:
    !###    COORDINATE_CONVERT_TO_RC_SP performs a coordinate transformation from a coordinate system identified by
    !###    COORDINATE_SYSTEM at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates for
    !###    single precision coordinates.
    !###  Parent-function: COORDINATE_CONVERT_TO_RC
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    REAL(SP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(SP) :: COORDINATE_CONVERT_TO_RC_SP(SIZE(X,1))
    !Local variables
    REAL(SP) :: FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_TO_RC_SP",ERR,ERROR,*999)
    
    COORDINATE_CONVERT_TO_RC_SP=0.0_SP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions",ERR,ERROR,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_TO_RC_SP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=X(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))
      CASE(3)
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_SP(3)=X(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN  
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(3)=X(1)*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
        COORDINATE_CONVERT_TO_RC_SP(1)=FOCUS*COSH(X(1))*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
        COORDINATE_CONVERT_TO_RC_SP(1)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(2)=FOCUS*SINH(X(1))*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_SP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("COORDINATE_CONVERT_TO_RC_SP")
    RETURN
999 CALL ERRORS("COORDINATE_CONVERT_TO_RC_SP",ERR,ERROR)
    CALL EXITS("COORDINATE_CONVERT_TO_RC_SP")
    RETURN 
  END FUNCTION COORDINATE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !

  !#### Generic-Function: COORDINATE_DELTA_CALCULATE
  !###  Description:
  !###    Calculates the difference (or delta) between two points in a coordinate system. Discontinuities for polar coordinate
  !###    systems are accounted for
  !###  Child-functions: COORDINATE_DELTA_CALCULATE_DP,COORDINATE_DELTA_CALCUALTE_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION COORDINATE_DELTA_CALCULATE_DP(COORDINATE_SYSTEM,X,Y,ERR,ERROR)
  
    !#### Function: COORDINATE_DELTA_CALCULATE_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    Calculates the difference (or detlta) between the point X and the point Y i.e., Y-X, in the given coordinate system.
    !###    0->2Pi discontinuities with polar coordinates are accounted for.
    !###  Parent-function: COORDINATE_DELTA_CALCULATE
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    REAL(DP), INTENT(IN) :: X(:)
    REAL(DP), INTENT(IN) :: Y(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: COORDINATE_DELTA_CALCULATE_DP(SIZE(X,1))
    !Local variables

    CALL ENTERS("COORDINATE_DELTA_CALCULATE_DP",ERR,ERROR,*999)

    COORDINATE_DELTA_CALCULATE_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions",ERR,ERROR,*999)

    IF(SIZE(X,1)/=SIZE(Y,1)) &
      & CALL FLAG_ERROR("Size of X is different to the size of Y",ERR,ERROR,*999)
   
    COORDINATE_DELTA_CALCULATE_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=Y(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)- &
      & X(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Do nothing
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("COORDINATE_DELTA_CALCULATE_DP")
    RETURN
999 CALL ERRORS("COORDINATE_DELTA_CALCULATE_DP",ERR,ERROR)
    CALL EXITS("COORDINATE_DELTA_CALCULATE_DP")
    RETURN 
  END FUNCTION COORDINATE_DELTA_CALCULATE_DP

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_METRICS_CALCULATE(COORDINATE_SYSTEM,JACOBIAN_TYPE,METRICS,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_METRICS_CALCUALTE
    !###  Description:
    !###    Calculates the covariant metric tensor GL(i,j), the contravariant metric tensor GU(i,J), the Jacobian and derivative
    !###    of the interpolated coordinate system (XI_i) with respect to the given coordinate (X_j) system (DXI_DX) at a point
    !###    (X - normally a Gauss point). 
    !###    Old-cmiss-name: XGMG

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: JACOBIAN_TYPE
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: METRICS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: mi,nj,ni,nu
    REAL(DP) :: DET_GL,DET_DX_DXI,FF,G1,G3,MU,R,RC,RCRC,RR
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("COORDINATE_METRICS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(METRICS)) THEN
        INTERPOLATED_POINT=>METRICS%INTERPOLATED_POINT
        IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
          IF(INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE>=FIRST_PART_DERIV) THEN
            
            !Calculate the derivatives of X with respect to XI
            DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
              METRICS%DX_DXI(1:METRICS%NUMBER_OF_X_DIMENSIONS,ni)=INTERPOLATED_POINT%VALUES(1:METRICS%NUMBER_OF_X_DIMENSIONS,nu)
            ENDDO !ni
            
            !Initialise the Jacobian and the metric tensors to the identity matrix
            METRICS%GL=0.0_DP
            METRICS%GU=0.0_DP
            DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
              METRICS%GL(ni,ni)=1.0_DP
              METRICS%GU(ni,ni)=1.0_DP
            ENDDO !i
            METRICS%JACOBIAN=0.0_DP
            
            !Calculate the covariant metric tensor GL(i,j)
            SELECT CASE(COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  DO nj=1,METRICS%NUMBER_OF_X_DIMENSIONS
                    METRICS%GL(mi,ni)=METRICS%GL(mi,ni)+METRICS%DX_DXI(nj,mi)*METRICS%DX_DXI(nj,ni)
                  ENDDO !nj
                ENDDO !ni
              ENDDO !mi
            CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
              R=INTERPOLATED_POINT%VALUES(1,1)
              RR=R*R
              IF(METRICS%NUMBER_OF_X_DIMENSIONS==2) THEN
                DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    METRICS%GL(mi,ni)=METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+RR*METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni)
                  ENDDO !ni
                ENDDO !mi
              ELSE IF(METRICS%NUMBER_OF_X_DIMENSIONS==3) THEN
                DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    METRICS%GL(mi,ni)=METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+RR*METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni)+ &
                      & METRICS%DX_DXI(3,mi)*METRICS%DX_DXI(3,ni)
                  ENDDO !ni
                ENDDO !mi
              ELSE
                LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(METRICS%NUMBER_OF_X_DIMENSIONS,"*",ERR,ERROR))// &
                  & " is an invalid number of dimensions for a cylindrical polar coordinate system"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
              R=INTERPOLATED_POINT%VALUES(1,1)
              RR=R*R
              RC=R*COS(INTERPOLATED_POINT%VALUES(3,1))
              RCRC=RC*RC          
              DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  METRICS%GL(mi,ni)=METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+RCRC*METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni)+ &
                    & RR*METRICS%DX_DXI(3,mi)*METRICS%DX_DXI(3,ni)
                ENDDO !ni
              ENDDO !mi
            CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
              IF(ABS(INTERPOLATED_POINT%VALUES(2,1))<ZERO_TOLERANCE) THEN
                CALL FLAG_WARNING("Mu is zero",ERR,ERROR,*999)
              ELSE
                FF=COORDINATE_SYSTEM%FOCUS*COORDINATE_SYSTEM%FOCUS
                R=INTERPOLATED_POINT%VALUES(1,1)
                MU=INTERPOLATED_POINT%VALUES(2,1)
                G1=FF*(SINH(R)*SINH(R)+SIN(MU)*SIN(MU))
                G3=FF*SINH(R)*SINH(R)*SIN(MU)*SIN(MU)
                IF(METRICS%NUMBER_OF_X_DIMENSIONS==2) THEN
                  DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                      METRICS%GL(mi,ni)=G1*(METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni))
                    ENDDO !ni
                  ENDDO !mi
                ELSE IF(METRICS%NUMBER_OF_X_DIMENSIONS==3) THEN
                  DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                      METRICS%GL(mi,ni)=G1*(METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni))+ &
                        & G3*METRICS%DX_DXI(3,mi)*METRICS%DX_DXI(3,ni)
                    ENDDO !ni
                  ENDDO !mi
                ELSE
                  LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(METRICS%NUMBER_OF_X_DIMENSIONS,"*",ERR,ERROR))// &
                    & " is an invalid number of dimensions for a prolate spheroidal coordinate system"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
              CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
                & " is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            
            !Calcualte the contravariant metric tensor
            CALL INVERT(METRICS%GL,METRICS%GU,DET_GL,ERR,ERROR,*999)
            
            !Calculate the derivatives of Xi with respect to X - DXI_DX
            IF(METRICS%NUMBER_OF_XI_DIMENSIONS==METRICS%NUMBER_OF_X_DIMENSIONS) THEN
              CALL INVERT(METRICS%DX_DXI,METRICS%DXI_DX,DET_DX_DXI,ERR,ERROR,*999)
            ELSE
              METRICS%DXI_DX=0.0_DP
            ENDIF
            
            !Calculate the Jacobian
            SELECT CASE(JACOBIAN_TYPE)
            CASE(COORDINATE_JACOBIAN_LINE_TYPE)
              METRICS%JACOBIAN=SQRT(ABS(METRICS%GL(1,1)))
              METRICS%JACOBIAN_TYPE=COORDINATE_JACOBIAN_LINE_TYPE
            CASE(COORDINATE_JACOBIAN_AREA_TYPE)
              IF(METRICS%NUMBER_OF_XI_DIMENSIONS==3) THEN
                METRICS%JACOBIAN=SQRT(ABS(DET_GL*METRICS%GU(3,3)))
              ELSE
                METRICS%JACOBIAN=SQRT(ABS(DET_GL))
              ENDIF
              METRICS%JACOBIAN_TYPE=COORDINATE_JACOBIAN_AREA_TYPE
            CASE(COORDINATE_JACOBIAN_VOLUME_TYPE)
              METRICS%JACOBIAN=SQRT(ABS(DET_GL))
              METRICS%JACOBIAN_TYPE=COORDINATE_JACOBIAN_VOLUME_TYPE
            CASE DEFAULT
              LOCAL_ERROR="The Jacobian type of "//TRIM(NUMBER_TO_VSTRING(JACOBIAN_TYPE,"*",ERR,ERROR))// &
                & " is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Metrics interpolated point has not been interpolated to include first derivatives",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Metrics interpolated point is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Metrics is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & COORDINATE_SYSTEM%TYPE)),ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",METRICS%NUMBER_OF_X_DIMENSIONS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",METRICS%NUMBER_OF_XI_DIMENSIONS,ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Location of metrics:",ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_X_DIMENSIONS,3,3,INTERPOLATED_POINT%VALUES(:,1), &
        & '("    X :",3(X,E13.6))','(7X,3(X,E13.6))',ERR,ERROR,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of X wrt Xi:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_X_DIMENSIONS,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS, &
        & 3,3,METRICS%DX_DXI,WRITE_STRING_MATRIX_NAME_AND_INDICIES,'("    dX_dXi','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Xi wrt X:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS,1,1,METRICS%NUMBER_OF_X_DIMENSIONS, &
        & 3,3,METRICS%DXI_DX,WRITE_STRING_MATRIX_NAME_AND_INDICIES,'("    dXi_dX','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Covariant metric tensor:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS, &
        & 3,3,METRICS%GL,WRITE_STRING_MATRIX_NAME_AND_INDICIES,'("    GL','(",I1,",:)',' :",3(X,E13.6))','(13X,3(X,E13.6))', &
        & ERR,ERROR,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Contravariant metric tensor:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS, &
        & 3,3,METRICS%GU,WRITE_STRING_MATRIX_NAME_AND_INDICIES,'("    GU','(",I1,",:)',' :",3(X,E13.6))','(13X,3(X,E13.6))', &
        & ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian type = ",METRICS%JACOBIAN_TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian = ",METRICS%JACOBIAN,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_METRICS_CALCULATE")
    RETURN
999 CALL ERRORS("COORDINATE_METRICS_CALCULATE",ERR,ERROR)
    CALL EXITS("COORDINATE_METRICS_CALCULATE")
    RETURN 1
  END SUBROUTINE COORDINATE_METRICS_CALCULATE

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_NORMAL_CALCULATE(COORDINATE_SYSTEM,REVERSE,X,N,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_NORMAL_CALCUALTE
    !###  Description:
    !###    Calculates the normal vector, N, at the point X. If REVERSE is true the reversed normal is returned 
    !###    Old-cmiss-name: NORMAL

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    LOGICAL, INTENT(IN) :: REVERSE
    REAL(DP), INTENT(IN) :: X(:,:)
    REAL(DP), INTENT(OUT) :: N(3)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_X_DIMENSIONS,d_s1,d_s2,d2_s1
    REAL(DP) :: LENGTH,R,TANGENT1(3),TANGENT2(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("COORDINATE_SYSTEM_NORMAL_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_X_DIMENSIONS=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        d_s1=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
        d_s2=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
        d2_s1=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(NUMBER_OF_X_DIMENSIONS==2) THEN
            TANGENT1(1)=X(1,d_s1)
            TANGENT1(2)=X(2,d_s1)
          ELSE IF(NUMBER_OF_X_DIMENSIONS==3) THEN
            TANGENT1(1)=X(1,d_s1)
            TANGENT1(2)=X(2,d_s1)
            TANGENT1(3)=X(3,d_s1)
            TANGENT2(1)=X(1,d_s2)
            TANGENT2(2)=X(2,d_s2)
            TANGENT2(3)=X(3,d_s2)
          ELSE
            LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(NUMBER_OF_X_DIMENSIONS,"*",ERR,ERROR))// &
              & " is an invalid number of dimensions to calculate a normal from in a rectangular cartesian coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
          R=X(1,1)
          IF(NUMBER_OF_X_DIMENSIONS==2) THEN
            TANGENT1(1)=X(1,d_s1)*COS(X(1,d_s1))-R*SIN(X(1,d_s1))*X(2,d_s1)
            TANGENT1(2)=X(2,d_s1)*SIN(X(1,d_s1))+R*COS(X(1,d_s1))*X(2,d_s1)
          ELSE IF(NUMBER_OF_X_DIMENSIONS==3) THEN
            TANGENT1(1)=X(1,d_s1)*COS(X(1,d_s1))-R*SIN(X(1,d_s1))*X(2,d_s1)
            TANGENT1(2)=X(2,d_s1)*SIN(X(1,d_s1))+R*COS(X(1,d_s1))*X(2,d_s1)
            TANGENT1(3)=X(3,d_s1)
            TANGENT2(1)=X(1,d_s2)*COS(X(1,d_s1))-R*SIN(X(1,d_s1))*X(2,d_s2)
            TANGENT2(2)=X(1,d_s2)*SIN(X(1,d_s1))+R*COS(X(1,d_s1))*X(2,d_s2)
            TANGENT2(3)=X(3,d_s2)
           ELSE
            LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(NUMBER_OF_X_DIMENSIONS,"*",ERR,ERROR))// &
              & " is an invalid number of dimensions to calculate a normal from in a rectangular cartesian coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF          
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          R=X(1,1)
          TANGENT1(1)=X(1,d_s1)*COS(X(1,d2_s1))*COS(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*COS(X(1,d_s1))*X(3,d_s1)- &
            &                 R*COS(X(1,d2_s1))*SIN(X(1,d_s1))*X(2,d_s1)
          TANGENT1(2)=X(1,d_s1)*COS(X(1,d2_s1))*SIN(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*SIN(X(1,d_s1))*X(3,d_s1)+ &
            &                 R*COS(X(1,d2_s1))*COS(X(1,d_s1))*X(2,d_s1)
          TANGENT1(3)=X(1,d_s1)*SIN(X(1,d2_s1))+R*COS(X(1,d2_s1))*X(3,d_s1)
          TANGENT2(1)=X(1,d_s2)*COS(X(1,d2_s1))*COS(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*COS(X(1,d_s1))*X(3,d_s2)- &
            &                 R*COS(X(1,d2_s1))*SIN(X(1,d_s1))*X(2,d_s2)
          TANGENT2(2)=X(1,d_s2)*COS(X(1,d2_s1))*SIN(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*SIN(X(1,d_s1))*X(3,d_s2)+ &
            &                 R*COS(X(1,d2_s1))*COS(X(1,d_s1))*X(2,d_s2)
          TANGENT2(3)=X(1,d_s2)*SIN(X(1,d2_s1))+R*COS(X(1,d2_s1))*X(3,d_s2)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
            & " is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(NUMBER_OF_X_DIMENSIONS==2) THEN
          N(1)=-TANGENT1(2)
          N(2)=TANGENT1(1)
          LENGTH=SQRT(N(1)*N(1)+N(2)*N(2))
          IF(ABS(LENGTH)<ZERO_TOLERANCE) CALL FLAG_ERROR("Zero normal vector length",ERR,ERROR,*999)
          IF(REVERSE) THEN
            N(1)=-N(1)/LENGTH
            N(2)=-N(2)/LENGTH
          ELSE            
            N(1)=N(1)/LENGTH
            N(2)=N(2)/LENGTH
          ENDIF
        ELSE
          N(1)=TANGENT1(2)*TANGENT2(3)-TANGENT1(3)*TANGENT2(2)
          N(2)=TANGENT1(3)*TANGENT2(1)-TANGENT1(1)*TANGENT2(3)
          N(3)=TANGENT1(1)*TANGENT2(2)-TANGENT1(2)*TANGENT2(1)
          LENGTH=SQRT(N(1)*N(1)+N(2)*N(2)+N(3)*N(3))
          IF(REVERSE) THEN
            N(1)=-N(1)/LENGTH
            N(2)=-N(2)/LENGTH
            N(3)=-N(3)/LENGTH
          ELSE            
            N(1)=N(1)/LENGTH
            N(2)=N(2)/LENGTH
            N(3)=N(3)/LENGTH
          ENDIF
        ENDIF        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & COORDINATE_SYSTEM%TYPE)),ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",NUMBER_OF_X_DIMENSIONS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_X_DIMENSIONS,3,3,X(:,1),'("  X         :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',ERR,ERROR,*999)      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_X_DIMENSIONS,3,3,TANGENT1,'("  Tangent 1 :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',ERR,ERROR,*999)
      IF(NUMBER_OF_X_DIMENSIONS==3) THEN
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_X_DIMENSIONS,3,3,TANGENT2,'("  Tangent 2 :",3(X,E13.6))', &
          & '(13X,3(X,E13.6))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_X_DIMENSIONS,3,3,N,'("  Normal    :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',ERR,ERROR,*999)            
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_NORMAL_CALCULATE")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_NORMAL_CALCULATE",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_NORMAL_CALCULATE")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_NORMAL_CALCULATE

  !
  !================================================================================================================================
  !

  PURE FUNCTION COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM)

    !#### Function: COORDINATE_SYSTEM_DIMENSION_GET
    !###  Type: INTEGER(INT
    !###  Description:
    !###    Gets the coordinate system dimension. Note: no error handling at the moment as there can be no errors for now.
    !###  See-Also: COORDINATE_SYSTEM_DIMENSION_SET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    !Function Result
    INTEGER(INTG) :: COORDINATE_SYSTEM_DIMENSION_GET
    !Local Variables

    COORDINATE_SYSTEM_DIMENSION_GET=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
    
    RETURN
  END FUNCTION COORDINATE_SYSTEM_DIMENSION_GET

  !
  !================================================================================================================================
  !

  FUNCTION COORDINATE_SYSTEM_FOCUS_GET(COORDINATE_SYSTEM,ERR,ERROR)

    !#### Function: COORDINATE_SYSTEM_FOCUS_GET
    !###  Type: REAL(DP)
    !###  Description:
    !###    Gets the coordinate system focus. 
    !###  See-Also: COORDINATE_SYSTEM_FOCUS_SET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function Result
    REAL(DP) :: COORDINATE_SYSTEM_FOCUS_GET
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_FOCUS_GET",ERR,ERROR,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      COORDINATE_SYSTEM_FOCUS_GET=COORDINATE_SYSTEM%FOCUS
    CASE DEFAULT
      CALL FLAG_ERROR("No focus defined for this coordinate system type",ERR,ERROR,*999)
    END SELECT
    
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_FOCUS_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_GET")
    RETURN
  END FUNCTION COORDINATE_SYSTEM_FOCUS_GET

  !
  !================================================================================================================================
  !

  PURE FUNCTION COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET(COORDINATE_SYSTEM)

    !#### Function: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Gets the coordinate system radial interpolation type. Note: no error handling at the moment as there can be no
    !###    errors for now.
    !###  See-Also: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    !Function Result
    INTEGER(INTG) :: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET
    !Local Variables

    COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET=COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE
    
    RETURN
  END FUNCTION COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET

  !
  !================================================================================================================================
  !

  PURE FUNCTION COORDINATE_SYSTEM_TYPE_GET(COORDINATE_SYSTEM)

    !#### Function: COORDINATE_SYSTEM_TYPE_GET
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Gets the coordinate system type. Note: no error handling at the moment as there can be no errors for now.
    !###  See-Also: COORDINATE_SYSTEM_TYPE_SET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    !Function Result
    INTEGER(INTG) :: COORDINATE_SYSTEM_TYPE_GET
    !Local Variables

    COORDINATE_SYSTEM_TYPE_GET=COORDINATE_SYSTEM%TYPE
    
    RETURN
  END FUNCTION COORDINATE_SYSTEM_TYPE_GET

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET_NUMBER(USER_NUMBER,DIMENSION,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_DIMENSION_SET_NUMBER
    !###  Description:
    !###    Sets/changes the dimension of the coordinate system identified by a USER_NUMBER.
    !###  See-Also: COORDINATE_SYSTEM_DIMENSION_GET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    INTEGER(INTG), INTENT(IN) :: DIMENSION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("COORDINATE_SYSTEM_DIMENSION_SET_NUMBER",ERR,ERROR,*999)

    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_DIMENSION_SET_PTR(COORDINATE_SYSTEM,DIMENSION,ERR,ERROR,*999)
    
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_SET_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DIMENSION_SET_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_SET_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET_PTR(COORDINATE_SYSTEM,DIMENSION,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_DIMENSION_SET_PTR
    !###  Description:
    !###    Sets/changes the dimension of the coordinate system identified by a ptr
    !###  See-Also: COORDINATE_SYSTEM_DIMENSION_GET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: DIMENSION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_DIMENSION_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(DIMENSION>=1.AND.DIMENSION<=3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
          IF(DIMENSION>=2.AND.DIMENSION<=3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_SET_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DIMENSION_SET_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_SET_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET_NUMBER(USER_NUMBER,FOCUS,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_FOCUS_SET_NUMBER
    !###  Description:
    !###    Sets/changes the focus of a coordinate system identified by a USER_NUMBER.
    !###  See-Also: COORDINATE_SYSTEM_FOCUS_GET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    REAL(DP), INTENT(IN) :: FOCUS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("COORDINATE_SYSTEM_FOCUS_SET_NUMBER",ERR,ERROR,*999)

    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_FOCUS_SET_PTR(COORDINATE_SYSTEM,FOCUS,ERR,ERROR,*999)
    
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_SET_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_FOCUS_SET_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_SET_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET_PTR(COORDINATE_SYSTEM,FOCUS,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_FOCUS_SET_PTR
    !###  Description:
    !###    Sets/changes the focus of a coordinate system identified by a pointer.
    !###  See-Also: COORDINATE_SYSTEM_FOCUS_GET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    REAL(DP), INTENT(IN) :: FOCUS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_FOCUS_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          IF(FOCUS>ZERO_TOLERANCE) THEN
            COORDINATE_SYSTEM%FOCUS=FOCUS
          ELSE
            CALL FLAG_ERROR("Focus is less than zero",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          IF(FOCUS>ZERO_TOLERANCE) THEN
            COORDINATE_SYSTEM%FOCUS=FOCUS
          ELSE
            CALL FLAG_ERROR("Focus is less than zero",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_SET_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_FOCUS_SET_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_SET_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER(USER_NUMBER,RADIAL_INTERPOLATION_TYPE,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER
    !###  Description:
    !###    Sets/changes the radial interpolation type of a coordinate system identified by a USER_NUMBER.
    !###  See-Also: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    INTEGER(INTG), INTENT(IN) :: RADIAL_INTERPOLATION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER",ERR,ERROR,*999)

    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR(COORDINATE_SYSTEM,RADIAL_INTERPOLATION_TYPE,ERR,ERROR,*999)
    
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR(COORDINATE_SYSTEM,RADIAL_INTERPOLATION_TYPE,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR
    !###  Description:
    !###    Sets/changes the radial interpolation type of a coordinate system identified by a pointer.
    !###  See-Also: COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: RADIAL_INTERPOLATION_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("COORDINATE_SYSTEM_INTERPOLATION_TYPE_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          SELECT CASE(RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a rectangular cartesian coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
          SELECT CASE(RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a cylindrical/spherical coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          SELECT CASE(RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a prolate spheroidal coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_TYPE_SET_NUMBER(USER_NUMBER,TYPE,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_TYPE_SET_NUMBER
    !###  Description:
    !###    Sets/changes the type of a coordinate system identified by a USER_NUMBER.
    !###  See-Also: COORDINATE_SYSTEM_TYPE_GET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    INTEGER(INTG), INTENT(IN) :: TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("COORDINATE_SYSTEM_TYPE_SET_NUMBER",ERR,ERROR,*999)

    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_TYPE_SET_PTR(COORDINATE_SYSTEM,TYPE,ERR,ERROR,*999)
    
    CALL EXITS("COORDINATE_SYSTEM_TYPE_SET_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_TYPE_SET_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_TYPE_SET_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_TYPE_SET_PTR(COORDINATE_SYSTEM,TYPE,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_TYPE_SET_PTR
    !###  Description:
    !###    Sets/changes the type of a coordinate system identified by a pointer.
    !###  See-Also: COORDINATE_SYSTEM_TYPE_GET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_TYPE_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
        CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_CYCLINDRICAL_POLAR_TYPE
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_SPHERICAL_POLAR_TYPE
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_PROLATE_SPHEROIDAL_TYPE
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_OBLATE_SPHEROIDAL_TYPE
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_TYPE_SET_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_TYPE_SET_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_TYPE_SET_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_TYPE_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET_NUMBER(USER_NUMBER,ORIGIN,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_ORIGIN_SET_NUMBER
    !###  Description:
    !###    Sets/changes the origin of a coordinate system identified by a number.
    !###  See-Also: COORDINATE_SYSTEM_TYPE_GET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    REAL(DP), INTENT(IN) :: ORIGIN(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("COORDINATE_SYSTEM_ORIGIN_SET_NUMBER",ERR,ERROR,*999)

    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_ORIGIN_SET_PTR(COORDINATE_SYSTEM,ORIGIN,ERR,ERROR,*999)
    
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_SET_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIGIN_SET_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_SET_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET_PTR(COORDINATE_SYSTEM,ORIGIN,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_ORIGIN_SET_PTR
    !###  Description:
    !###    Sets/changes the origin of a coordinate system identified by a pointer.
    !###  See-Also: COORDINATE_SYSTEM_ORIGIN_GET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    REAL(DP), INTENT(IN) :: ORIGIN(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_ORIGIN_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        IF(SIZE(ORIGIN)==3) THEN
          COORDINATE_SYSTEM%ORIGIN=ORIGIN
        ELSE
          CALL FLAG_ERROR("The origin must have exactly 3 components",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_SET_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIGIN_SET_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_SET_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER(USER_NUMBER,ORIENTATION,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER
    !###  Description:
    !###    Sets/changes the orientation of a coordinate system identified by a USER_NUMBER.
    !###  See-Also: COORDINATE_SYSTEM_ORIENTATION_GET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    REAL(DP), INTENT(IN) :: ORIENTATION(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER",ERR,ERROR,*999)

    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_ORIENTATION_SET_PTR(COORDINATE_SYSTEM,ORIENTATION,ERR,ERROR,*999)
    
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET_PTR(COORDINATE_SYSTEM,ORIENTATION,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_ORIENTATION_SET_PTR
    !###  Description:
    !###    Sets/changes the orientation of a coordinate system identified by a pointer.
    !###  See-Also: COORDINATE_SYSTEM_ORIENTATION_GET

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    REAL(DP), INTENT(IN) :: ORIENTATION(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_ORIENTATION_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished",ERR,ERROR,*999)
      ELSE
        IF(SIZE(ORIENTATION,1)==3.AND.SIZE(ORIENTATION,2)==3) THEN
          !!TODO: Check orientation matrix vectors are orthogonal to each other etc.
          COORDINATE_SYSTEM%ORIENTATION=ORIENTATION
        ELSE
          CALL FLAG_ERROR("The orientation matrix must have exactly 3x3 components",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_SET_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIENTATION_SET_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_SET_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_CREATE_START(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_CREATE_START
    !###  Description:
    !###    Starts the creation of and initialises a new coordinate system. The coordinate system is 3D rectangular cartesian by
    !###    default. Defaults may be changed by COORDINATE_SYSTEM_SET_xxx calls.
    !###  See-Also: COORDINATE_SYSTEM_DESTROY,COORDINATE_SYSTEM_CREATE_FINISH

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: NEW_COORDINATE_SYSTEM
    TYPE(COORDINATE_SYSTEM_PTR_TYPE), POINTER :: NEW_COORDINATE_SYSTEMS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_COORDINATE_SYSTEM)
    NULLIFY(NEW_COORDINATE_SYSTEMS)

    CALL ENTERS("COORDINATE_SYSTEM_CREATE_START",ERR,ERROR,*999)

    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      LOCAL_ERROR="Coordinate system number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
        & " has already been created"
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ELSE
      ALLOCATE(NEW_COORDINATE_SYSTEM,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordinate system",ERR,ERROR,*999)
      
      NEW_COORDINATE_SYSTEM%USER_NUMBER=USER_NUMBER
      NEW_COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED=.FALSE.
      NEW_COORDINATE_SYSTEM%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      NEW_COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
      NEW_COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=3
      NEW_COORDINATE_SYSTEM%FOCUS=1.0_DP    
      NEW_COORDINATE_SYSTEM%ORIGIN=(/0.0_DP,0.0_DP,0.0_DP/)
      NEW_COORDINATE_SYSTEM%ORIENTATION=RESHAPE(&
        & (/1.0_DP,0.0_DP,0.0_DP, &
        &   0.0_DP,1.0_DP,0.0_DP, &
        &   0.0_DP,0.0_DP,1.0_DP/), &
        & (/3,3/))
      
      ALLOCATE(NEW_COORDINATE_SYSTEMS(COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS+1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordinate systems",ERR,ERROR,*999)
      DO coord_system_idx=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
        NEW_COORDINATE_SYSTEMS(coord_system_idx)%PTR=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR
      ENDDO !coord_system_idx
      NEW_COORDINATE_SYSTEMS(COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS+1)%PTR=>NEW_COORDINATE_SYSTEM
      DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
      COORDINATE_SYSTEMS%COORDINATE_SYSTEMS=>NEW_COORDINATE_SYSTEMS
      COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS+1
      
      COORDINATE_SYSTEM=>NEW_COORDINATE_SYSTEM

    ENDIF
        
    CALL EXITS("COORDINATE_SYSTEM_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_COORDINATE_SYSTEM)) DEALLOCATE(NEW_COORDINATE_SYSTEM)
    IF(ASSOCIATED(NEW_COORDINATE_SYSTEMS)) DEALLOCATE(NEW_COORDINATE_SYSTEMS)
    NULLIFY(COORDINATE_SYSTEM)
    CALL ERRORS("COORDINATE_SYSTEM_CREATE_START",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_CREATE_START")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_CREATE_FINISH
    !###  Description:
    !###    Finishes the creation of a new coordinate system.

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: coord_system_idx

    CALL ENTERS("COORDINATE_SYSTEM_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of coordinate systems = ", &
        & COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS,ERR,ERROR,*999)
      DO coord_system_idx=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system : ",coord_system_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number = ", &
          & COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Type = ", &
          & COORDINATE_SYSTEM_TYPE_STRING(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR%TYPE),ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ", &
          & COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
      ENDDO !coord_system_idx
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_CREATE_FINISH")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_DESTROY_NUMBER(USER_NUMBER,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_DESTROY_NUMBER
    !###  Description:
    !###    Destroys a coordinate system with the the given USER_NUMBER.

    !Argument variables
    INTEGER(INTG) :: USER_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: coord_system_idx,new_coord_system_idx
    LOGICAL :: FOUND
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(COORDINATE_SYSTEM_PTR_TYPE), POINTER :: NEW_COORDINATE_SYSTEMS(:)

    CALL ENTERS("COORDINATE_SYSTEM_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(USER_NUMBER==0) THEN
      CALL FLAG_ERROR("Cannot destroy the global coordinate system",ERR,ERROR,*999)
    ELSE
      FOUND=.FALSE.
      new_coord_system_idx=0
      ALLOCATE(NEW_COORDINATE_SYSTEMS(COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS-1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordianate systems",ERR,ERROR,*999)
      DO coord_system_idx=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
        IF(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          FOUND=.TRUE.
          COORDINATE_SYSTEM=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR
        ELSE
          new_coord_system_idx=new_coord_system_idx+1
          NEW_COORDINATE_SYSTEMS(new_coord_system_idx)%PTR=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR
        ENDIF
      ENDDO !coord_system_idx
      IF(FOUND) THEN
        DEALLOCATE(COORDINATE_SYSTEM)
        DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
        COORDINATE_SYSTEMS%COORDINATE_SYSTEMS=>NEW_COORDINATE_SYSTEMS
        COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS-1
      ELSE
        DEALLOCATE(NEW_COORDINATE_SYSTEMS)
        CALL FLAG_ERROR("Coordinate system number to destroy does not exist",ERR,ERROR,*999)
      ENDIF
    ENDIF
        
    CALL EXITS("COORDINATE_SYSTEM_DESTROY_NUMBER")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DESTROY_NUMBER")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_DESTROY_PTR(COORDINATE_SYSTEM,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_DESTROY_PTR
    !###  Description:
    !###    Destroys a coordinate system given by a pointer to the coordinate system

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: coord_system_no,new_coord_system_no
    LOGICAL :: FOUND
    TYPE(COORDINATE_SYSTEM_PTR_TYPE), POINTER :: NEW_COORDINATE_SYSTEMS(:)

    CALL ENTERS("COORDINATE_SYSTEM_DESTROY_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%USER_NUMBER==0) THEN
        CALL FLAG_ERROR("Cannot destroy the global coordinate system",ERR,ERROR,*999)
      ELSE
        FOUND=.FALSE.
        new_coord_system_no=0
        ALLOCATE(NEW_COORDINATE_SYSTEMS(COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS-1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordianate systems",ERR,ERROR,*999)
        DO coord_system_no=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
          IF(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_no)%PTR%USER_NUMBER==COORDINATE_SYSTEM%USER_NUMBER) THEN
            FOUND=.TRUE.
          ELSE
            new_coord_system_no=new_coord_system_no+1
            NEW_COORDINATE_SYSTEMS(new_coord_system_no)%PTR=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_no)%PTR
          ENDIF
        ENDDO !coord_system_no
        IF(FOUND) THEN
          DEALLOCATE(COORDINATE_SYSTEM)
          DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
          COORDINATE_SYSTEMS%COORDINATE_SYSTEMS=>NEW_COORDINATE_SYSTEMS
          COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS-1
        ELSE
          DEALLOCATE(NEW_COORDINATE_SYSTEMS)
          CALL FLAG_ERROR("Coordinate system number to destroy does not exist",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("COORDINATE_SYSTEM_DESTROY_PTR")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DESTROY_PTR",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DESTROY_PTR")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DESTROY_PTR

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: DXZ
  !###  Description:
  !###    Calculates DX(:)/DZ(I) at X, where Z(I) are rectangular 
  !###    Cartesian and X(:) are curvilinear coordinates defined by 
  !###    COORDINATE_SYSTEM.
  !###  Child-functions: DXZ_DP,DXZ_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION DXZ_DP(COORDINATE_SYSTEM,I,X,ERR,ERROR)
  
    !#### Function: DXZ_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    Calculates DX(:)/DZ(I) at X, where Z(I) are rectangular 
    !###    Cartesian and X(:) are curvilinear coordinates defined by 
    !###    COORDINATE_SYSTEM for double precision coordinates.
    !###  Parent-function: DXZ
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: I
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: DXZ_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: RD,FOCUS

    CALL ENTERS("DXZ_DP",ERR,ERROR,*999)

    DXZ_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions", &
      & ERR,ERROR,*999)
   
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      IF(I>0.AND.I<=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
        DXZ_DP(I)=1.0_DP
      ELSE
        CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=COS(X(2))
          DXZ_DP(2)=-SIN(X(2))/X(1)
        CASE(2)
          DXZ_DP(1)=SIN(X(2))
          DXZ_DP(2)=COS(X(2))/X(1)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      CASE(3)
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=COS(X(2))
          DXZ_DP(2)=-SIN(X(2))/X(1)
          DXZ_DP(3)=0.0_DP
        CASE(2)
          DXZ_DP(1)=SIN(X(2))
          DXZ_DP(2)=COS(X(2))/X(1)
          DXZ_DP(3)=0.0_DP
        CASE(3)
          DXZ_DP(1)=0.0_DP
          DXZ_DP(2)=0.0_DP
          DXZ_DP(3)=1.0_DP
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=COS(X(2))*COS(X(3))
          DXZ_DP(2)=-SIN(X(2))/(X(1)*COS(X(3)))
          DXZ_DP(3)=-COS(X(2))*SIN(X(3))/X(1)
        CASE(2)
          DXZ_DP(1)=SIN(X(2))*COS(X(3))
          DXZ_DP(2)=COS(X(2))/(X(1)*COS(X(3)))
          DXZ_DP(3)=-SIN(X(2))*SIN(X(3))/X(1)
        CASE(3)
          DXZ_DP(1)=SIN(X(3))
          DXZ_DP(2)=0.0_DP
          DXZ_DP(3)=COS(X(3))/X(1)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        RD=FOCUS*(COSH(X(1))*COSH(X(1))-COS(X(2))*COS(X(2)))
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=SINH(X(1))*COS(X(2))/RD
          DXZ_DP(2)=-COSH(X(1))*SIN(X(2))/RD
          DXZ_DP(3)=0.0_DP
        CASE(2)
          DXZ_DP(1)=COSH(X(1))*SIN(X(2))*COS(X(3))/RD
          DXZ_DP(2)=SINH(X(1))*COS(X(2))*COS(X(3))/RD
          DXZ_DP(3)=-SIN(X(3))/(FOCUS*SINH(X(1))*SIN(X(2)))
        CASE(3)
          DXZ_DP(1)=COSH(X(1))*SIN(X(2))*SIN(X(3))/RD
          DXZ_DP(2)=SINH(X(1))*COS(X(2))*SIN(X(3))/RD
          DXZ_DP(3)=COS(X(3))/(FOCUS*SINH(X(1))*SIN(X(2)))
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("DXZ_DP")
    RETURN
999 CALL ERRORS("DXZ_DP",ERR,ERROR)
    CALL EXITS("DXZ_DP")
    RETURN 
  END FUNCTION DXZ_DP

  !
  !================================================================================================================================
  !

  !#### Generic-Function: D2ZX
  !###  Description:
  !###    Calculates D2Z(:)/DX(I)DX(J) at X(:), where Z(:) are rectangalar
  !###    Cartesian and X(I) and X(J) are curvilinear coordinates defined by 
  !###    COORDINATE_SYSTEM.
  !###  Child-functions: D2ZX_SP,D2ZX_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION D2ZX_DP(COORDINATE_SYSTEM,I,J,X,ERR,ERROR)
  
    !#### Function: D2ZX_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    Calculates D2Z(:)/DX(I)DX(J) at X(:), where Z(:) are rectangalar
    !###    Cartesian and X(I) and X(J) are curvilinear coordinates defined by 
    !###    COORDINATE_SYSTEM.
    !###  Parent-function: D2ZX
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: I,J
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: D2ZX_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS

    CALL ENTERS("D2ZX_DP",ERR,ERROR,*999)

    D2ZX_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions", &
      & ERR,ERROR,*999)
   
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      D2ZX_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=0.0_DP
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
          CASE(2)
            D2ZX_DP(1)=-X(1)*COS(X(2))
            D2ZX_DP(2)=-X(1)*SIN(X(2))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      CASE(3)
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-X(1)*COS(X(2))
            D2ZX_DP(2)=-X(1)*SIN(X(2))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(3)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-SIN(X(2))*COS(X(3))
            D2ZX_DP(2)=COS(X(2))*COS(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=-COS(X(2))*SIN(X(3))
            D2ZX_DP(2)=-SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=COS(X(3))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-SIN(X(2))*COS(X(3))
            D2ZX_DP(2)=COS(X(2))*COS(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-X(1)*COS(X(2))*COS(X(3))
            D2ZX_DP(2)=-X(1)*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=X(1)*SIN(X(2))*SIN(X(3))
            D2ZX_DP(2)=-X(1)*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(3)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-COS(X(2))*SIN(X(3))
            D2ZX_DP(2)=-SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=COS(X(3))
          CASE(2)
            D2ZX_DP(1)=X(1)*SIN(X(2))*SIN(X(3))
            D2ZX_DP(2)=-X(1)*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=-X(1)*COS(X(2))*COS(X(3))
            D2ZX_DP(2)=-X(1)*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=-X(1)*SIN(X(3))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=FOCUS*COSH(X(1))*COS(X(2))          
            D2ZX_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          CASE(2)
            D2ZX_DP(1)=-FOCUS*SINH(X(1))*SIN(X(2))
            D2ZX_DP(2)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-FOCUS*SINH(X(1))*SIN(X(2))
            D2ZX_DP(2)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
          CASE(2)
            D2ZX_DP(1)=-FOCUS*COSH(X(1))*COS(X(2))
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=-FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE(3)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          CASE(2)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=-FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid j value",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("D2ZX_DP")
    RETURN
999 CALL ERRORS("D2ZX_DP",ERR,ERROR)
    CALL EXITS("D2ZX_DP")
    RETURN 
  END FUNCTION D2ZX_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: DZX
  !###  Description:
  !###    Calculates DZ(:)/DX(J) at X, where Z(:) are rectangalar
  !###    Cartesian and X(J) are curvilinear coordinates defined by 
  !###    COORDINATE_SYSTEM.
  !###  Child-functions: DZX_DP,DZX_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION DZX_DP(COORDINATE_SYSTEM,I,X,ERR,ERROR)
  
    !#### Function: DZX_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    Calculates DZ(:)/DX(I) at X, where Z(:) are rectangalar
    !###    Cartesian and X(I) are curvilinear coordinates defined by 
    !###    COORDINATE_SYSTEM.
    !###  Parent-function: DZX
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: I
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: DZX_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS

    CALL ENTERS("DZX_DP",ERR,ERROR,*999)

    DZX_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions", &
      & ERR,ERROR,*999)
   
   SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      IF(I>0.AND.I<=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
        DZX_DP(I)=1.0_DP
      ELSE
        CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=COS(X(2))
          DZX_DP(2)=SIN(X(2))
        CASE(2)
          DZX_DP(1)=-X(1)*SIN(X(2))
          DZX_DP(2)=X(1)*COS(X(2))
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      CASE(3)
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=COS(X(2))
          DZX_DP(2)=SIN(X(2))
          DZX_DP(3)=0.0_DP
        CASE(2)
          DZX_DP(1)=-X(1)*SIN(X(2))
          DZX_DP(2)=X(1)*COS(X(2))
          DZX_DP(3)=0.0_DP
        CASE(3)
          DZX_DP(1)=0.0_DP
          DZX_DP(2)=0.0_DP
          DZX_DP(3)=1.0_DP
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=COS(X(2))*COS(X(3))
          DZX_DP(2)=SIN(X(2))*COS(X(3))
          DZX_DP(3)=SIN(X(3))
        CASE(2)
          DZX_DP(1)=-X(1)*SIN(X(2))*COS(X(3))
          DZX_DP(2)=X(1)*COS(X(2))*COS(X(3))
          DZX_DP(3)=0.0_DP
        CASE(3)
          DZX_DP(1)=-X(1)*COS(X(2))*SIN(X(3))
          DZX_DP(2)=-X(1)*SIN(X(2))*SIN(X(3))
          DZX_DP(3)=X(1)*COS(X(3))
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=FOCUS*SINH(X(1))*COS(X(2))
          DZX_DP(2)=FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          DZX_DP(3)=FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
        CASE(2)
          DZX_DP(1)=-FOCUS*COSH(X(1))*SIN(X(2))
          DZX_DP(2)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          DZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
        CASE(3)
          DZX_DP(1)=0.0_DP
          DZX_DP(2)=-FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          DZX_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          DZX_DP(2)=FOCUS*COSH(X(1))*SIN(X(2))
          DZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
        CASE(2)
          DZX_DP(1)=-FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          DZX_DP(2)=FOCUS*SINH(X(1))*COS(X(2))
          DZX_DP(3)=-FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
        CASE(3)
          DZX_DP(1)=-FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
          DZX_DP(2)=0.0_DP
          DZX_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid i value",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
      ENDIF
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("DZX_DP")
    RETURN
999 CALL ERRORS("DZX_DP",ERR,ERROR)
    CALL EXITS("DZX_DP")
    RETURN 
  END FUNCTION DZX_DP

  !
  !================================================================================================================================
  !

  !!TODO:: CHANGE THIS TO A FUNCTION
  
  !#### Generic-subroutine: COORDINATE_DERIVATIVE_CONVERT_TO_RC
  !###  Description:
  !###    COORDINATE_DERIVATIVE_CONVERT_TO_RC performs a coordinate transformation from a
  !###    coordinate system identified by COORDINATE at the point X
  !###    with coordinates/derivatives X(nj,nu) to the point with 
  !###    coordinates/derivatives Z(nj) in rectangular cartesian coordinates
  !###  Child-subroutines: COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP,COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP(COORDINATE_SYSTEM,PART_DERIV_TYPE,X,Z,&
    & ERR,ERROR,*)
  
    !#### Subroutine: COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP
    !###  Description:
    !###    COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP performs a coordinate transformation from a
    !###    coordinate system identified by COORDINATE_SYSTEM at the point X
    !###    with coordinates/derivatives X(nj,nu) to the point with 
    !###    coordinates/derivatives Z(nj) in rectangular cartesian coordinates
    !###    for double precision coordinates.
    !###  Parent-function: COORDINATE_DERIVATIVE_CONVERT_TO_RC
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_TYPE
    REAL(DP), INTENT(IN) :: X(:,:)
    REAL(DP), INTENT(OUT) :: Z(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    REAL(DP) :: FOCUS
    
    CALL ENTERS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP",ERR,ERROR,*999)

!!TODO: change all second index X(:,?) numbers to their apropriate constant
!!as defined in constants e.g. X(1,2) == X(1,PART_DERIV_S1)
    
    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions", &
      & ERR,ERROR,*999)
    
    IF(SIZE(X,1)==SIZE(Z,1)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        IF(SIZE(X,2)>=PART_DERIV_TYPE) THEN
          Z=X(:,PART_DERIV_TYPE)
        ELSE
          CALL FLAG_ERROR("Invalid derivative type",ERR,ERROR,*999)
        ENDIF
      CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
        SELECT CASE(PART_DERIV_TYPE)
        CASE(NO_PART_DERIV)
          Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
          IF(ERR/=0) GOTO 999
        CASE(PART_DERIV_S1)
          IF(SIZE(X,2)>=2) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
              Z(3)=X(3,2) !d(z)/d(s1)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S1_S1)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PART_DERIV_S2)
          IF(SIZE(X,2)>=4) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
              Z(3)=X(3,4) !d(z)/d(s2)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S2_S2)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PART_DERIV_S1_S2)
          IF(SIZE(X,2)>=6) THEN
            SELECT CASE(SIZE(X,1))
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
              Z(3)=X(3,6) !d2(z)/d(s1)d(s2)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S3)
          IF(SIZE(X,2)>=7) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,7)-X(1,1)*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=SIN(X(2,1))*X(1,7)+X(1,1)*COS(X(2,1))*X(2,7) !d(y)/d(s3)
              Z(3)=X(3,7) !d(z)/d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S3_S3)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PART_DERIV_S1_S3)
          IF(SIZE(X,2)>=9) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,9)-X(1,2)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,9)) !d2(x)/d(s1)d(s3)
              Z(2)=SIN(X(2,1))*X(1,9)+X(1,2)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,9)) !d2(y)/d(s1)d(s3)
              Z(3)=X(3,9) !d2(z)/d(s1)d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S2_S3)
          IF(SIZE(X,2)>=10) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,10)-X(1,4)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,4)+X(1,1)*COS(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,10)) !d2(x)/d(s2)d(s3)
              Z(2)=SIN(X(2,1))*X(1,10)+X(1,4)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,4)-X(1,1)*SIN(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,10)) !d2(y)/d(s2)d(s3)
              Z(3)=X(3,10) !d2(z)/d(s2)d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S1_S2_S3)
          IF(SIZE(X,2)>=11) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)  
              Z(1)=-COS(X(2,1))*X(2,2)*X(2,4)*X(1,7)-&
                & SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))-&
                & SIN(X(2,1))*X(2,2)*X(1,10)+COS(X(2,1))*X(1,11)-&
                & COS(X(2,1))*X(2,2)*X(1,4)*X(2,7)-& 
                & SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COS(X(2,1))*X(1,2)+SIN(X(2,1))*X(1,1)*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*COS(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(-SIN(X(2,1))*X(1,2)-X(1,1)*& 
                & COS(X(2,1))*X(2,2))*X(2,10)-X(1,1)*&
                & SIN(X(2,1))*X(2,11) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=-SIN(X(2,1))*X(2,2)*X(2,4)*X(1,7)+&
                & COS(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & COS(X(2,1))*X(2,2)*X(1,10)+SIN(X(2,1))*X(1,11)-&
                & SIN(X(2,1))*X(2,2)*X(1,4)*X(2,7)+&
                & COS(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SIN(X(2,1))*X(1,2)-X(1,1)*COS(X(2,1))*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*SIN(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(COS(X(2,1))*X(1,2)-X(1,1)*&
                & SIN(X(2,1))*X(2,2))*X(2,10)+X(1,1)*&
                & COS(X(2,1))*X(2,11) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=X(3,11) !d3(z)/d(s1)d(s2)d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
        END SELECT
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
            IF(ERR/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
        ENDIF
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
            IF(ERR/=0) GOTO 999
          CASE(PART_DERIV_S1)
            IF(SIZE(X,2)>=2) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,2) !d(y)/d(s1)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,2) !d(z)/d(s1)
            ELSE
              CALL FLAG_ERROR("Invalid derivative type",ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S1_S1)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE(PART_DERIV_S2)
            IF(SIZE(X,2)>=4) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,4)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4) !d(y)/d(s2)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,4)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4) !d(z)/d(s2)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S2_S2)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE(PART_DERIV_S1_S2)
            IF(SIZE(X,2)>=6) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,6)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,4))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,6)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,4)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,4))) !d2(x)/d(s1)d(s2)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,6)+X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,4)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,4))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,6)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))) !d2(y)/d(s1)d(s2)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,6)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,6)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))) !d2(z)/d(s1)d(s2)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S3)
            IF(SIZE(X,2)>=7) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7) !d(y)/d(s3)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7) !d(z)/d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S3_S3)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE(PART_DERIV_S1_S3)
            IF(SIZE(X,2)>=9) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,9)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,9)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s1)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,9)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))* X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s1)d(s3) 
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,9)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s1)d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S2_S3)
            IF(SIZE(X,2)>=10) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,10)+&
                & X(1,4)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,10)-&
                & X(2,4)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s2)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,10)-&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s2)d(s3)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,10)+&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s2)d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S1_S2_S3)
            IF(SIZE(X,2)>=11) THEN
              Z(1)=FOCUS*((SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,4)*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*(X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,4)*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,10)+&
                & SINH(X(1,1))*COS(X(2,1))*X(1,11)+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(1,4)*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SINH(X(1,1))*COS(X(2,1))*X(1,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(2,4)*X(2,7)-&
                & COSH(X(1,1))*COS(X(2,1))*(X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,10)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,11)) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,11))+&
                & FOCUS*((-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(3,10)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,11)) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,11))+&
                & FOCUS*((SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(3,10)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,11)) !d3(z)/d(s1)d(s2)d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
          END SELECT
        ENDIF
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
            IF(ERR/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("The sizes of the vectors X and Z do not match",&
        & ERR,ERROR,*999)
    ENDIF

    CALL EXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP")
    RETURN
999 CALL ERRORS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP",ERR,ERROR)
    CALL EXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP")
    RETURN 1
  END SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP(COORDINATE_SYSTEM,PART_DERIV_TYPE,X,Z,&
    & ERR,ERROR,*)
  
    !#### Subroutine: COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP
    !###  Description:
    !###    COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP performs a coordinate transformation from a
    !###    coordinate system identified by COORDINATE_SYSTEM at the point X
    !###    with coordinates/derivatives X(nj,nu) to the point with 
    !###    coordinates/derivatives Z(nj) in rectangular cartesian coordinates
    !###    for single precision coordinates.
    !###  Parent-function: COORDINATE_DERIVATIVE_CONVERT_TO_RC
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_TYPE
    REAL(SP), INTENT(IN) :: X(:,:)
    REAL(SP), INTENT(OUT) :: Z(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    REAL(SP) :: FOCUS
    
    CALL ENTERS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP",ERR,ERROR,*999)
    
!!TODO: change all second index X(:,?) numbers to their apropriate constant
!!as defined in constants e.g. X(1,2) == X(1,PART_DERIV_S1)

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions", &
      & ERR,ERROR,*999)
    
    IF(SIZE(X,1)==SIZE(Z,1)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        IF(SIZE(X,2)>=PART_DERIV_TYPE) THEN
          Z=X(:,PART_DERIV_TYPE)
        ELSE
          CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
        ENDIF
      CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
        SELECT CASE(PART_DERIV_TYPE)
        CASE(NO_PART_DERIV)
          Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
          IF(ERR/=0) GOTO 999
        CASE(PART_DERIV_S1)
          IF(SIZE(X,2)>=2) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
              Z(3)=X(3,2) !d(z)/d(s1)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S1_S1)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PART_DERIV_S2)
          IF(SIZE(X,2)>=4) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
              Z(3)=X(3,4) !d(z)/d(s2)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S2_S2)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PART_DERIV_S1_S2)
          IF(SIZE(X,2)>=6) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
              Z(3)=X(3,6) !d2(z)/d(s1)d(s2)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S3)
          IF(SIZE(X,2)>=7) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,7)-X(1,1)*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=SIN(X(2,1))*X(1,7)+X(1,1)*COS(X(2,1))*X(2,7) !d(y)/d(s3)
              Z(3)=X(3,7) !d(z)/d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S3_S3)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(PART_DERIV_S1_S3)
          IF(SIZE(X,2)>=9) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,9)-X(1,2)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,9)) !d2(x)/d(s1)d(s3)
              Z(2)=SIN(X(2,1))*X(1,9)+X(1,2)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,9)) !d2(y)/d(s1)d(s3)
              Z(3)=X(3,9) !d2(z)/d(s1)d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S2_S3)
          IF(SIZE(X,2)>=10) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,10)-X(1,4)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,4)+X(1,1)*COS(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,10)) !d2(x)/d(s2)d(s3)
              Z(2)=SIN(X(2,1))*X(1,10)+X(1,4)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,4)-X(1,1)*SIN(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,10)) !d2(y)/d(s2)d(s3)
              Z(3)=X(3,10) !d2(z)/d(s2)d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE(PART_DERIV_S1_S2_S3)
          IF(SIZE(X,2)>=11) THEN
            SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)  
              Z(1)=-COS(X(2,1))*X(2,2)*X(2,4)*X(1,7)-&
                & SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))-&
                & SIN(X(2,1))*X(2,2)*X(1,10)+COS(X(2,1))*X(1,11)-&
                & COS(X(2,1))*X(2,2)*X(1,4)*X(2,7)-& 
                & SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COS(X(2,1))*X(1,2)+SIN(X(2,1))*X(1,1)*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*COS(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(-SIN(X(2,1))*X(1,2)-X(1,1)*& 
                & COS(X(2,1))*X(2,2))*X(2,10)-X(1,1)*&
                & SIN(X(2,1))*X(2,11) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=-SIN(X(2,1))*X(2,2)*X(2,4)*X(1,7)+&
                & COS(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & COS(X(2,1))*X(2,2)*X(1,10)+SIN(X(2,1))*X(1,11)-&
                & SIN(X(2,1))*X(2,2)*X(1,4)*X(2,7)+&
                & COS(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SIN(X(2,1))*X(1,2)-X(1,1)*COS(X(2,1))*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*SIN(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(COS(X(2,1))*X(1,2)-X(1,1)*&
                & SIN(X(2,1))*X(2,2))*X(2,10)+X(1,1)*&
                & COS(X(2,1))*X(2,11) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=X(3,11) !d3(z)/d(s1)d(s2)d(s3)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Not enough X derivatives supplied",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
        END SELECT
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
            IF(ERR/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
        ENDIF
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
            IF(ERR/=0) GOTO 999
          CASE(PART_DERIV_S1)
            IF(SIZE(X,2)>=2) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,2) !d(y)/d(s1)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,2) !d(z)/d(s1)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied", &
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S1_S1)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE(PART_DERIV_S2)
            IF(SIZE(X,2)>=4) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,4)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4) !d(y)/d(s2)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,4)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4) !d(z)/d(s2)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S2_S2)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE(PART_DERIV_S1_S2)
            IF(SIZE(X,2)>=6) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,6)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,4))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,6)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,4)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,4))) !d2(x)/d(s1)d(s2)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,6)+X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,4)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,4))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,6)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))) !d2(y)/d(s1)d(s2)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,6)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,6)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))) !d2(z)/d(s1)d(s2)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S3)
            IF(SIZE(X,2)>=7) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7) !d(y)/d(s3)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7) !d(z)/d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S3_S3)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE(PART_DERIV_S1_S3)
            IF(SIZE(X,2)>=9) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,9)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,9)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s1)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,9)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))* X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s1)d(s3) 
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,9)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s1)d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S2_S3)
            IF(SIZE(X,2)>=10) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,10)+&
                & X(1,4)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,10)-&
                & X(2,4)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s2)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,10)-&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s2)d(s3)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,10)+&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s2)d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE(PART_DERIV_S1_S2_S3)
            IF(SIZE(X,2)>=11) THEN
              Z(1)=FOCUS*((SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,4)*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*(X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,4)*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,10)+&
                & SINH(X(1,1))*COS(X(2,1))*X(1,11)+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(1,4)*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SINH(X(1,1))*COS(X(2,1))*X(1,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(2,4)*X(2,7)-&
                & COSH(X(1,1))*COS(X(2,1))*(X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,10)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,11)) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,11))+&
                & FOCUS*((-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(3,10)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,11)) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,11))+&
                & FOCUS*((SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(3,10)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,11)) !d3(z)/d(s1)d(s2)d(s3)
            ELSE
              CALL FLAG_ERROR("Not enough X derivatives supplied",&
                & ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
          END SELECT
        ENDIF
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),ERR,ERROR)
            IF(ERR/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid partial derivative type",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Invalid number of coordinates",ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid coordinate type",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("The sizes of the vectors X and Z do not match",&
        & ERR,ERROR,*999)
    ENDIF

    CALL EXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP")
    RETURN
999 CALL ERRORS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP",ERR,ERROR)
    CALL EXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP")
    RETURN 1
  END SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_DERIVATIVE_NORM(COORDINATE_SYSTEM,PART_DERIV_INDEX,INTERPOLATED_POINT,DERIV_NORM,ERR,ERROR,*)
  
    !#### Subroutine: COORDINATE_DERIVATIVE_NORM
    !###  Description:
    !###    Calculates the norm of a derivative in a coordinate system identified by COORDINATE_SYSTEM at the given interpolated
    !###    point and returns the value in NORM for single precision coordinates. PART_DERIV_INDEX is used to select the
    !###    appropriate partial derivative (i.e., wrt S1, S2 or S3) to normalise.
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_INDEX
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: DERIV_NORM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: component_idx,NUMBER_OF_COMPONENTS
    REAL(DP) :: FOCUS,SL,SM
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("COORDINATE_DERIVATIVE_NORM",ERR,ERROR,*999)

    DERIV_NORM=0.0_DP
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
        IF(INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE>=FIRST_PART_DERIV) THEN
          IF(PART_DERIV_INDEX<=INTERPOLATED_POINT%MAX_PARTIAL_DERIVATIVE_INDEX) THEN
            NUMBER_OF_COMPONENTS=SIZE(INTERPOLATED_POINT%VALUES,1)
            SELECT CASE(PART_DERIV_INDEX)
            CASE(PART_DERIV_S1,PART_DERIV_S2,PART_DERIV_S3)
              SELECT CASE(COORDINATE_SYSTEM%TYPE)
              CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  DERIV_NORM=DERIV_NORM+INTERPOLATED_POINT%VALUES(component_idx,PART_DERIV_INDEX)**2
                ENDDO !component_index
              CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
                IF(NUMBER_OF_COMPONENTS==2) THEN
                  DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                    & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2
                ELSE IF(NUMBER_OF_COMPONENTS==3) THEN
                  DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                    & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2+INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX)**2
                ELSE
                  LOCAL_ERROR="The number of components for the interpolated point of "// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                    & " is invalid for a cylindrical polar coordinate system"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
                DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                  & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX)*COS(INTERPOLATED_POINT%VALUES(3,1)))**2+ &
                  & (INTERPOLATED_POINT%VALUES(1,1)*INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX))**2
              CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
                FOCUS=COORDINATE_SYSTEM%FOCUS
                SL=SINH(INTERPOLATED_POINT%VALUES(1,1))
                SM=SIN(INTERPOLATED_POINT%VALUES(2,1))
                DERIV_NORM=FOCUS*FOCUS*((SL*SL+SM*SM)*(INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+ &
                  & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2)+(SL*SM*INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX))**2
              CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
                CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
                  & " is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              DERIV_NORM=SQRT(DERIV_NORM)
            CASE DEFAULT
              LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PART_DERIV_INDEX,"*",ERR,ERROR))// &
                & " is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PART_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The interpolated point has a maximum number of partial derivatives of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATED_POINT%MAX_PARTIAL_DERIVATIVE_INDEX,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The point has not been interpolated to include first derivative values",ERR,ERROR,*999)
        ENDIF          
      ELSE
        CALL FLAG_ERROR("Interpolated point is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_DERIVATIVE_NORM")
    RETURN
999 CALL ERRORS("COORDINATE_DERIVATIVE_NORM",ERR,ERROR)
    CALL EXITS("COORDINATE_DERIVATIVE_NORM")
    RETURN 1
  END SUBROUTINE COORDINATE_DERIVATIVE_NORM

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,PARTIAL_DERIVATIVE_INDEX,VALUE,ERR,ERROR,*)
  
    !#### Subroutine: COORDINATE_INTERPOLATION_ADJUST
    !###  Description:
    !###    Adjusts the interpolation for non-rectangular cartesian coordinate systems.
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX
    REAL(DP), INTENT(INOUT) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    REAL(DP) :: COSHX,CSS,D,DES,FOCUS,R,SS,SINHX,THETA
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("COORDINATE_INTERPOLATION_ADJUST",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        !Do nothing
      CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
        SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
        CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
          !Do nothing
        CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          R=SQRT(VALUE)
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=R
          ELSE
            VALUE=VALUE/(2.0_DP*R)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
            & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a cylindrical coordinate system"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
        CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
          !Do nothing
        CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          R=VALUE**(1.0_DP/3.0_DP)
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=R
          ELSE
            VALUE=VALUE/(3.0_DP*R*R)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
            & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a cylindrical/spherical coordinate system"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
        CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
          !Do nothing
        CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SS=VALUE/(FOCUS*FOCUS)
          SINHX=SQRT(SS)
          COSHX=SQRT(1.0_DP+SS)
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=LOG(SINHX+COSHX)
          ELSE
            VALUE=VALUE/(2.0_DP*FOCUS*FOCUS*SINHX*COSHX)
          ENDIF
        CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
          FOCUS=COORDINATE_SYSTEM%FOCUS
          CSS=VALUE/(FOCUS**3.0_DP)
          DES=CSS*CSS-4.0_DP/27.0_DP
          IF(DES>0.0_DP) THEN
            D=((CSS+SQRT(DES))/2.0_DP)**(1.0_DP/3.0_DP)
            COSHX=D+1.0_DP/(3.0_DP*D)
          ELSE
            THETA=ACOS(CSS*SQRT(27.0_DP)/2.0_DP)
            COSHX=2.0_DP/SQRT(3.0_DP)*COS(THETA/3.0_DP)
          ENDIF
          SINHX=SQRT(ABS(COSHX*COSHX-1.0_DP))
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=LOG(SINHX+COSHX)
          ELSE
            VALUE=VALUE/((3.0_DP*COSHX*COSHX-1.0_DP)*SINHX)/(FOCUS**3.0_DP)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
            & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a prolate spheroidal coordinate system"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
          & " is invalid"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("COORDINATE_INTERPOLATION_ADJUST")
    RETURN
999 CALL ERRORS("COORDINATE_INTERPOLATION_ADJUST",ERR,ERROR)
    CALL EXITS("COORDINATE_INTERPOLATION_ADJUST")
    RETURN 1
  END SUBROUTINE COORDINATE_INTERPOLATION_ADJUST

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_INTERPOLATION_PARAMETERS_ADJUST(COORDINATE_SYSTEM,INTERPOLATION_PARAMETERS,ERR,ERROR,*)
  
    !#### Subroutine: COORDINATE_INTERPOLATION_PARAMETERS_ADJUST
    !###  Description:
    !###    Adjusts the interpolation parameters for non-rectangular cartesian coordinate systems.
    
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST",ERR,ERROR,*999)

!!TODO: Tidy up element parameters for non-rc coordinate systems. See bottom of XPXE and ZPZE.
    
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          !Do nothing
        CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a cylindrical coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a spherical coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a prolate spheroidal coordinate system"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
            & " is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Interpolation parameters is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST")
    RETURN
999 CALL ERRORS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST",ERR,ERROR)
    CALL EXITS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST")
    RETURN 1
  END SUBROUTINE COORDINATE_INTERPOLATION_PARAMETERS_ADJUST

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEM_USER_NUMBER_FIND
    !###  Description:
    !###   Returns a pointer to the coordinate system identified by USER_NUMBER. If a coordinate system with that number is not
    !###   found then COORDINATE_SYSTEM is set to NULL.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    
    CALL ENTERS("COORDINATE_SYSTEM_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(USER_NUMBER==0) THEN
      COORDINATE_SYSTEM=>GLOBAL_COORDINATE_SYSTEM
    ELSE
      NULLIFY(COORDINATE_SYSTEM)
      coord_system_idx=1
      DO WHILE(coord_system_idx<=COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS.AND..NOT.ASSOCIATED(COORDINATE_SYSTEM))
        IF(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          COORDINATE_SYSTEM=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR
        ELSE
          coord_system_idx=coord_system_idx+1
        ENDIF
      ENDDO
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEMS_FINALISE(ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEMS_FINALISE
    !###  Description:
    !###   Finalises the coordinate systems and destroys all coordinate systems.
    !###  See-Also: COORDINATE_SYSTEMS_INITIALISE

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    
    CALL ENTERS("COORDINATE_SYSTEMS_FINALISE",ERR,ERROR,*999)

    DO coord_system_idx=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
      DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR)
    ENDDO !coord_system_idx
    DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
    COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=0
    NULLIFY(GLOBAL_COORDINATE_SYSTEM)
    
    CALL EXITS("COORDINATE_SYSTEMS_FINALISE")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEMS_FINALISE",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEMS_FINALISE")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEMS_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE COORDINATE_SYSTEMS_INITIALISE(ERR,ERROR,*)

    !#### Subroutine: COORDINATE_SYSTEMS_INITIALISE
    !###  Description:
    !###   Initialises the coordinate systems and creates the world coordinate system.
    !###  See-Also: COORDINATE_SYSTEMS_FINALISE

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("COORDINATE_SYSTEMS_INITIALISE",ERR,ERROR,*999)

    !Create the default RC World cooordinate system
    ALLOCATE(GLOBAL_COORDINATE_SYSTEM,STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global coordinate system",ERR,ERROR,*999)
    GLOBAL_COORDINATE_SYSTEM%USER_NUMBER=0
    GLOBAL_COORDINATE_SYSTEM%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    GLOBAL_COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=3
    GLOBAL_COORDINATE_SYSTEM%FOCUS=1.0_DP    
    GLOBAL_COORDINATE_SYSTEM%ORIGIN=(/0.0_DP,0.0_DP,0.0_DP/)
    GLOBAL_COORDINATE_SYSTEM%ORIENTATION=RESHAPE(&
      & (/1.0_DP,0.0_DP,0.0_DP, &
      &   0.0_DP,1.0_DP,0.0_DP, &
      &   0.0_DP,0.0_DP,1.0_DP/), &
      & (/3,3/))    
    !Store the global system in the list of coordinate systems
    ALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate coordinate systems",ERR,ERROR,*999)
    COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=1
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR=>GLOBAL_COORDINATE_SYSTEM
    
    CALL EXITS("COORDINATE_SYSTEMS_INITIALISE")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEMS_INITIALISE",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEMS_INITIALISE")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEMS_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE COORDINATE_ROUTINES
