!> \file
!> $Id$
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
!> Contributor(s): Kumar Mithraratne
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
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
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
  INTEGER(INTG), PARAMETER :: COORDINATE_CYLINDRICAL_POLAR_TYPE=2 !<Cylindrical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
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
  
  !Interfaces

  !>COORDINATE_CONVERT_FROM_RC performs a coordinate transformation from a rectangular cartesian coordinate at the point with
  !>coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by COORDINATE_SYSTEM.
  INTERFACE COORDINATE_CONVERT_FROM_RC
    MODULE PROCEDURE COORDINATE_CONVERT_FROM_RC_DP
    MODULE PROCEDURE COORDINATE_CONVERT_FROM_RC_SP
  END INTERFACE !COORDINATE_CONVERT_FROM_RC

  !>COORDINATE_CONVERT_TO_RC performs a coordinate transformation from a coordinate system identified by COORDINATE_SYSTEM
  !>at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates.
  INTERFACE COORDINATE_CONVERT_TO_RC
    MODULE PROCEDURE COORDINATE_CONVERT_TO_RC_DP
    MODULE PROCEDURE COORDINATE_CONVERT_TO_RC_SP
  END INTERFACE !COORDINATE_CONVERT_TO_RC

  !>Calculates the difference (or delta) between two points in a coordinate system. Discontinuities for polar coordinate
  !>systems are accounted for
  INTERFACE COORDINATE_DELTA_CALCULATE
    MODULE PROCEDURE COORDINATE_DELTA_CALCULATE_DP
    !MODULE PROCEDURE COORDINATE_DELTA_CALCULATE_SP
  END INTERFACE !COORDINATE_DELTA_CALCULATE

  !>Calculates DX(:)/DZ(I) at X, where Z(I) are rectangular cartesian and X(:) are curvilinear coordinates defined by COORDINATE_SYSTEM. \todo CHANGE NAME TO SOMETHING MORE MEANINGFULL?
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

  PUBLIC COORDINATE_RECTANGULAR_CARTESIAN_TYPE,COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE, &
    & COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE

  PUBLIC COORDINATE_NO_RADIAL_INTERPOLATION_TYPE,COORDINATE_RADIAL_INTERPOLATION_TYPE, &
    & COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE,COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE
  
  PUBLIC COORDINATE_JACOBIAN_LINE_TYPE,COORDINATE_JACOBIAN_AREA_TYPE,COORDINATE_JACOBIAN_VOLUME_TYPE
  
  PUBLIC COORDINATE_SYSTEM_TYPE_STRING

  PUBLIC COORDINATE_CONVERT_FROM_RC,COORDINATE_CONVERT_TO_RC,COORDINATE_DELTA_CALCULATE,COORDINATE_DERIVATIVE_NORM, &
    & COORDINATE_INTERPOLATION_ADJUST,COORDINATE_INTERPOLATION_PARAMETERS_ADJUST, &
    & COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE,COORDINATE_METRICS_CALCULATE
  
  PUBLIC COORDINATE_SYSTEM_DIMENSION_GET,COORDINATE_SYSTEM_DIMENSION_SET

  PUBLIC COORDINATE_SYSTEM_FOCUS_GET,COORDINATE_SYSTEM_FOCUS_SET

  PUBLIC COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET,COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET
  
  PUBLIC COORDINATE_SYSTEM_TYPE_GET,COORDINATE_SYSTEM_TYPE_SET

  PUBLIC COORDINATE_SYSTEM_ORIGIN_GET,COORDINATE_SYSTEM_ORIGIN_SET
  
  PUBLIC COORDINATE_SYSTEM_ORIENTATION_GET,COORDINATE_SYSTEM_ORIENTATION_SET

  PUBLIC COORDINATE_SYSTEM_CREATE_START,COORDINATE_SYSTEM_CREATE_FINISH

  PUBLIC COORDINATE_SYSTEM_DESTROY

  PUBLIC COORDINATE_DERIVATIVE_CONVERT_TO_RC

  PUBLIC COORDINATE_SYSTEM_USER_NUMBER_FIND
  
  PUBLIC COORDINATE_SYSTEMS_INITIALISE,COORDINATE_SYSTEMS_FINALISE
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Performs a coordinate transformation from a rectangular cartesian coordinate at the point with coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by COORDINATE_SYSTEM for double precision coordinates.
  FUNCTION COORDINATE_CONVERT_FROM_RC_DP(COORDINATE_SYSTEM,Z,ERR,ERROR)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to perform the conversion on
    REAL(DP), INTENT(IN) :: Z(:) !<The rectangular cartesian coordiantes to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: COORDINATE_CONVERT_FROM_RC_DP(SIZE(Z,1))
    !Local variables
    REAL(DP) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_FROM_RC_DP",ERR,ERROR,*999)

    COORDINATE_CONVERT_FROM_RC_DP=0.0_DP

    IF(SIZE(Z,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of Z is less than the number of dimensions.",ERR,ERROR,*999)
    
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_FROM_RC_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=Z(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(1),Z(2))
      CASE(3)
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(1),Z(2))
        COORDINATE_CONVERT_FROM_RC_DP(3)=Z(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
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
          CALL FLAG_WARNING("Put A8=0 since ABS(A8)>1.",ERR,ERROR,*999)
        ENDIF
        IF((Z(3)==0.0_DP).OR.(A6==0.0_DP).OR.(A7==0.0_DP)) THEN
          A9=0.0_DP
        ELSE
          IF(ABS(A6*A7)>0.0_DP) THEN
            A9=Z(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0_DP
            CALL FLAG_WARNING("Put A9=0 since A6*A7=0.",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type.",ERR,ERROR,*999)
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

  !>COORDINATE_CONVERT_FROM_RC_SP performs a coordinate transformation from a rectangular cartesian coordinate at the
  !>point with coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by
  !>COORDINATE_SYSTEM for single precision coordinates.
  FUNCTION COORDINATE_CONVERT_FROM_RC_SP(COORDINATE_SYSTEM,Z,ERR,ERROR)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to convert from RC to
    REAL(SP), INTENT(IN) :: Z(:) !<The coordinate to convert from rectangular cartesian
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(SP) :: COORDINATE_CONVERT_FROM_RC_SP(SIZE(Z,1))
    !Local variables
    REAL(SP) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_FROM_RC_SP",ERR,ERROR,*999)

    COORDINATE_CONVERT_FROM_RC_SP=0.0_SP
    
    IF(SIZE(Z,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of Z is less than the number of dimensions.",ERR,ERROR,*999)
    
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_FROM_RC_SP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=Z(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(1),Z(2))
      CASE(3)
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(1),Z(2))
        COORDINATE_CONVERT_FROM_RC_SP(3)=Z(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
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
          CALL FLAG_WARNING("Put A8=0 since ABS(A8)>1.",ERR,ERROR,*999)
        ENDIF
        IF((Z(3)==0.0_SP).OR.(A6==0.0_SP).OR.(A7==0.0_SP)) THEN
          A9=0.0_SP
        ELSE
          IF(ABS(A6*A7)>0.0_SP) THEN
            A9=Z(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0_SP
            CALL FLAG_WARNING("Put A9=0 since A6*A7=0.",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type.",ERR,ERROR,*999)
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

  !>COORDINATE_CONVERT_TO_RC_DP performs a coordinate transformation from a coordinate system identified by
  !>COORDINATE_SYSTEM at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates for
  !>double precision coordinates.
  FUNCTION COORDINATE_CONVERT_TO_RC_DP(COORDINATE_SYSTEM,X,ERR,ERROR)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to convert to rectangular cartesian
    REAL(DP), INTENT(IN) :: X(:) !<The coordiante to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error coode
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: COORDINATE_CONVERT_TO_RC_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_TO_RC_DP",ERR,ERROR,*999)
    
    COORDINATE_CONVERT_TO_RC_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions.",ERR,ERROR,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_TO_RC_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=X(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))
      CASE(3)
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_DP(3)=X(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN  
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(3)=X(1)*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        COORDINATE_CONVERT_TO_RC_DP(1)=FOCUS*COSH(X(1))*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        COORDINATE_CONVERT_TO_RC_DP(1)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type.",ERR,ERROR,*999)
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

  !>COORDINATE_CONVERT_TO_RC_SP performs a coordinate transformation from a coordinate system identified by
  !>COORDINATE_SYSTEM at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates for
  !>single precision coordinates.
  FUNCTION COORDINATE_CONVERT_TO_RC_SP(COORDINATE_SYSTEM,X,ERR,ERROR)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to convert to rectangular cartesian
    REAL(SP), INTENT(IN) :: X(:) !<The coordinate to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(SP) :: COORDINATE_CONVERT_TO_RC_SP(SIZE(X,1))
    !Local variables
    REAL(SP) :: FOCUS
    
    CALL ENTERS("COORDINATE_CONVERT_TO_RC_SP",ERR,ERROR,*999)
    
    COORDINATE_CONVERT_TO_RC_SP=0.0_SP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions.",ERR,ERROR,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_TO_RC_SP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=X(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
      CASE(2)
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))
      CASE(3)
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_SP(3)=X(3)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN  
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(3)=X(1)*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
        COORDINATE_CONVERT_TO_RC_SP(1)=FOCUS*COSH(X(1))*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS==3) THEN
        FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
        COORDINATE_CONVERT_TO_RC_SP(1)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(2)=FOCUS*SINH(X(1))*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_SP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
      ELSE
        CALL FLAG_ERROR("Invalid number of coordinates.",ERR,ERROR,*999)
      ENDIF
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type.",ERR,ERROR,*999)
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

  !>Calculates the difference (or detlta) between the point X and the point Y i.e., Y-X, in the given coordinate system.
  !>0->2Pi discontinuities with polar coordinates are accounted for.
  FUNCTION COORDINATE_DELTA_CALCULATE_DP(COORDINATE_SYSTEM,X,Y,ERR,ERROR)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to calculate the delta for
    REAL(DP), INTENT(IN) :: X(:) !<The first coordinate
    REAL(DP), INTENT(IN) :: Y(:) !<The second coordinate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: COORDINATE_DELTA_CALCULATE_DP(SIZE(X,1))
    !Local variables

    CALL ENTERS("COORDINATE_DELTA_CALCULATE_DP",ERR,ERROR,*999)

    COORDINATE_DELTA_CALCULATE_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) &
      & CALL FLAG_ERROR("Size of X is less than the number of dimensions.",ERR,ERROR,*999)

    IF(SIZE(X,1)/=SIZE(Y,1)) &
      & CALL FLAG_ERROR("Size of X is different to the size of Y.",ERR,ERROR,*999)
   
    COORDINATE_DELTA_CALCULATE_DP(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)=Y(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)- &
      & X(1:COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS)
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Do nothing
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid coordinate type.",ERR,ERROR,*999)
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

  !>Calculates the covariant metric tensor GL(i,j), the contravariant metric tensor GU(i,J), the Jacobian and derivative of the interpolated coordinate system (XI_i) with respect to the given coordinate (X_j) system (DXI_DX) at a point (X - normally a Gauss point). Old cmiss name: XGMG
  SUBROUTINE COORDINATE_METRICS_CALCULATE(COORDINATE_SYSTEM,JACOBIAN_TYPE,METRICS,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to calculate the metrics for
    INTEGER(INTG), INTENT(IN) :: JACOBIAN_TYPE !<The type of Jacobian to calculate \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: METRICS !<A pointer to the metrics to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mi,ni,nu
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
              SELECT CASE(METRICS%NUMBER_OF_X_DIMENSIONS)
              CASE(1)
                DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    METRICS%GL(mi,ni)=METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)                    
                  ENDDO !ni
                ENDDO !mi
              CASE(2)
                DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    METRICS%GL(mi,ni)=METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+ &
                      & METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni)
                  ENDDO !ni
                ENDDO !mi
              CASE(3)
                DO mi=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                  DO ni=1,METRICS%NUMBER_OF_XI_DIMENSIONS
                    METRICS%GL(mi,ni)=METRICS%DX_DXI(1,mi)*METRICS%DX_DXI(1,ni)+ &
                      & METRICS%DX_DXI(2,mi)*METRICS%DX_DXI(2,ni)+ &
                      & METRICS%DX_DXI(3,mi)*METRICS%DX_DXI(3,ni)
                  ENDDO !ni
                ENDDO !mi
              CASE DEFAULT
                LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(METRICS%NUMBER_OF_X_DIMENSIONS,"*",ERR,ERROR))// &
                  & " is an invalid number of dimensions for a rectangular cartesian coordinate system."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
                  & " is an invalid number of dimensions for a cylindrical polar coordinate system."
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
                CALL FLAG_WARNING("Mu is zero.",ERR,ERROR,*999)
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
                    & " is an invalid number of dimensions for a prolate spheroidal coordinate system."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            
            !Calcualte the contravariant metric tensor
            CALL INVERT(METRICS%GL(1:METRICS%NUMBER_OF_XI_DIMENSIONS,1:METRICS%NUMBER_OF_XI_DIMENSIONS), &
              & METRICS%GU(1:METRICS%NUMBER_OF_XI_DIMENSIONS,1:METRICS%NUMBER_OF_XI_DIMENSIONS),DET_GL, &
              & ERR,ERROR,*999)
            
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
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Metrics interpolated point has not been interpolated to include first derivatives.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Metrics interpolated point is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Metrics is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & COORDINATE_SYSTEM%TYPE)),ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",METRICS%NUMBER_OF_X_DIMENSIONS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",METRICS%NUMBER_OF_XI_DIMENSIONS,ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Location of metrics:",ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_X_DIMENSIONS,3,3,INTERPOLATED_POINT%VALUES(:,1), &
        & '("    X           :",3(X,E13.6))','(17X,3(X,E13.6))',ERR,ERROR,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of X wrt Xi:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_X_DIMENSIONS,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS, &
        & 3,3,METRICS%DX_DXI,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dX_dXi','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Xi wrt X:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS,1,1,METRICS%NUMBER_OF_X_DIMENSIONS, &
        & 3,3,METRICS%DXI_DX,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dXi_dX','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Covariant metric tensor:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS, &
        & 3,3,METRICS%GL,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    GL','(",I1,",:)','     :",3(X,E13.6))','(17X,3(X,E13.6))', &
        & ERR,ERROR,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Contravariant metric tensor:",ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS,1,1,METRICS%NUMBER_OF_XI_DIMENSIONS, &
        & 3,3,METRICS%GU,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    GU','(",I1,",:)','     :",3(X,E13.6))','(17X,3(X,E13.6))', &
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

  !>Calculates the normal vector, N, at the point X. IF REVERSE is true the reversed normal is returned. Old-cmiss-name: NORMAL
  SUBROUTINE COORDINATE_SYSTEM_NORMAL_CALCULATE(COORDINATE_SYSTEM,REVERSE,X,N,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to calculate the normal for
    LOGICAL, INTENT(IN) :: REVERSE !<If .TRUE. the reversed normal is returned.
    REAL(DP), INTENT(IN) :: X(:,:) !<The coordinate and it's derivatives to calcualte the normal at
    REAL(DP), INTENT(OUT) :: N(3) !<On exit, the normal vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_X_DIMENSIONS,d_s1,d_s2,d2_s1
    REAL(DP) :: LENGTH,R,TANGENT1(3),TANGENT2(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("COORDINATE_SYSTEM_NORMAL_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
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
              & " is an invalid number of dimensions to calculate a normal from in a rectangular cartesian coordinate system."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
              & " is an invalid number of dimensions to calculate a normal from in a rectangular cartesian coordinate system."
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
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(NUMBER_OF_X_DIMENSIONS==2) THEN
          N(1)=-TANGENT1(2)
          N(2)=TANGENT1(1)
          LENGTH=SQRT(N(1)*N(1)+N(2)*N(2))
          IF(ABS(LENGTH)<ZERO_TOLERANCE) CALL FLAG_ERROR("Zero normal vector length.",ERR,ERROR,*999)
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
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
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

  !>Gets the coordinate system dimension. 
  SUBROUTINE COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,NUMBER_OF_DIMENSIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the dimension for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_DIMENSIONS !<On return, the number of dimensions in the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_DIMENSION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        NUMBER_OF_DIMENSIONS=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
      ELSE
        CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
   
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DIMENSION_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_GET")
    RETURN 1

  END SUBROUTINE COORDINATE_SYSTEM_DIMENSION_GET

  !
  !================================================================================================================================
  !

  !>Finalises a coordinate system and deallocates all memory. 
  SUBROUTINE COORDINATE_SYSTEM_FINALISE(COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      DEALLOCATE(COORDINATE_SYSTEM)
    ENDIF
   
    CALL EXITS("COORDINATE_SYSTEM_FINALISE")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_FINALISE",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_FINALISE")
    RETURN 1

  END SUBROUTINE COORDINATE_SYSTEM_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system focus. 
  SUBROUTINE COORDINATE_SYSTEM_FOCUS_GET(COORDINATE_SYSTEM,FOCUS,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the focus for
    REAL(DP), INTENT(OUT) :: FOCUS !<On return, the focus of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_FOCUS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          FOCUS=COORDINATE_SYSTEM%FOCUS
        CASE DEFAULT
          CALL FLAG_ERROR("No focus defined for this coordinate system type.",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_FOCUS_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_GET")
    RETURN 1
    
  END SUBROUTINE COORDINATE_SYSTEM_FOCUS_GET

  !
  !================================================================================================================================
  !


  !>Gets the coordinate system radial interpolation type. 
  SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET(COORDINATE_SYSTEM,RADIAL_INTERP_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to get the radial interpolation for
    INTEGER(INTG), INTENT(OUT) :: RADIAL_INTERP_TYPE !<On return, the radial interpolation type for the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
          RADIAL_INTERP_TYPE=COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE
        CASE DEFAULT
          CALL FLAG_ERROR("No radial interpolation type defined for this coordinate system interpolation.",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET")
    RETURN 1
    
  END SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system type. 
  SUBROUTINE COORDINATE_SYSTEM_TYPE_GET(COORDINATE_SYSTEM,SYSTEM_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the type for
    INTEGER(INTG), INTENT(OUT) :: SYSTEM_TYPE !<On return, the type for the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        SYSTEM_TYPE=COORDINATE_SYSTEM%TYPE
      ELSE
        CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_TYPE_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_TYPE_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_TYPE_GET")
    RETURN 1

  END SUBROUTINE COORDINATE_SYSTEM_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension of the coordinate system. \see OPENCMISS::CMISSCoordinateSystemDimensionSet
  SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,DIMENSION,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer the coordinate system to set the dimension for
    INTEGER(INTG), INTENT(IN) :: DIMENSION !<The dimension to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_DIMENSION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(DIMENSION>=1.AND.DIMENSION<=3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          IF(DIMENSION>=2.AND.DIMENSION<=3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS=DIMENSION
          ELSE
            CALL FLAG_ERROR("Invalid number of dimensions.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type.",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_SET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DIMENSION_SET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DIMENSION_SET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the focus of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemFocusSet
  SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET(COORDINATE_SYSTEM,FOCUS,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to set the focus for
    REAL(DP), INTENT(IN) :: FOCUS !<The focus to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_FOCUS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          IF(FOCUS>ZERO_TOLERANCE) THEN
            COORDINATE_SYSTEM%FOCUS=FOCUS
          ELSE
            CALL FLAG_ERROR("Focus is less than zero.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          IF(FOCUS>ZERO_TOLERANCE) THEN
            COORDINATE_SYSTEM%FOCUS=FOCUS
          ELSE
            CALL FLAG_ERROR("Focus is less than zero.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type.",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_SET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_FOCUS_SET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_FOCUS_SET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the radial interpolation type of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemRadialInterpolationTypeSet
  SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET(COORDINATE_SYSTEM,RADIAL_INTERPOLATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set the interpolation type for
    INTEGER(INTG), INTENT(IN) :: RADIAL_INTERPOLATION_TYPE !<The interpolation type to set \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("COORDINATE_SYSTEM_INTERPOLATION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          SELECT CASE(RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a rectangular cartesian coordinate system."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
          SELECT CASE(RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
            COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE=COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a cylindrical/spherical coordinate system."
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
              & " is invalid for a prolate spheroidal coordinate system."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type.",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_RADIAL_INTERPOLATION_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemTypeSet
  SUBROUTINE COORDINATE_SYSTEM_TYPE_SET(COORDINATE_SYSTEM,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The coordinate system type to set \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_CYLINDRICAL_POLAR_TYPE
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_SPHERICAL_POLAR_TYPE
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_PROLATE_SPHEROIDAL_TYPE
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_OBLATE_SPHEROIDAL_TYPE
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid coordinate system type.",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_TYPE_SET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_TYPE_SET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_TYPE_SET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Returns the origin of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOriginGet
  SUBROUTINE COORDINATE_SYSTEM_ORIGIN_GET(COORDINATE_SYSTEM,ORIGIN,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the origin for
    REAL(DP), INTENT(OUT) :: ORIGIN(:) !<On return, the origin of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_ORIGIN_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        IF(SIZE(ORIGIN)>=3) THEN
          ORIGIN(1:3)=COORDINATE_SYSTEM%ORIGIN
        ELSE
          CALL FLAG_ERROR("The origin must have >= 3 components.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIGIN_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_GET")
    RETURN 1
    
  END SUBROUTINE COORDINATE_SYSTEM_ORIGIN_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOriginSet
  SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET(COORDINATE_SYSTEM,ORIGIN,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to set the origin for
    REAL(DP), INTENT(IN) :: ORIGIN(:) !<The origin to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_ORIGIN_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(ORIGIN)==3) THEN
          COORDINATE_SYSTEM%ORIGIN=ORIGIN
        ELSE
          CALL FLAG_ERROR("The origin must have exactly 3 components.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_SET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIGIN_SET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIGIN_SET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET

  !
  !================================================================================================================================
  !

  !>Returns the orientation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOrientationSets/changesGet
  SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_GET(COORDINATE_SYSTEM,ORIENTATION,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the orientation for
    REAL(DP), INTENT(OUT) :: ORIENTATION(:,:) !<On return, the orientation of the coordinate system
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_ORIENTATION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        IF(SIZE(ORIENTATION,1)>=3.AND.SIZE(ORIENTATION,2)>=3) THEN
          ORIENTATION(1:3,1:3)=COORDINATE_SYSTEM%ORIENTATION
        ELSE
          CALL FLAG_ERROR("The orientation matrix must have >= 3x3 components.",ERR,ERROR,*999)
        ENDIF
      ELSE
         CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_GET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIENTATION_GET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_GET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the orientation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOrientationSet
  SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET(COORDINATE_SYSTEM,ORIENTATION,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to set the orientation for
    REAL(DP), INTENT(IN) :: ORIENTATION(:,:) !<The orientation to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COORDINATE_SYSTEM_ORIENTATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
        CALL FLAG_ERROR("Coordinate system has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(ORIENTATION,1)==3.AND.SIZE(ORIENTATION,2)==3) THEN
!!TODO: \todo Check orientation matrix vectors are orthogonal to each other etc.
          COORDINATE_SYSTEM%ORIENTATION=ORIENTATION
        ELSE
          CALL FLAG_ERROR("The orientation matrix must have exactly 3x3 components.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_SET")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_ORIENTATION_SET",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_ORIENTATION_SET")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET
  
  !
  !================================================================================================================================
  !

  !>Starts the creation of and initialises a new coordinate system. \see OPENCMISS::CMISSCoordinateSystemCreateStart
  !>The default values of the COORDINATE_SYSTEM's attributes are:
  !>- TYPE: 1 (COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
  !>- RADIAL_INTERPOLATION_TYPE: 0 (COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
  !>- Dimensions: 3
  !>- Focus: 1.0
  !>- Origin: (0.0,0.0,0.0)
  !>- Oritention: ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
  SUBROUTINE COORDINATE_SYSTEM_CREATE_START(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number for the created coordinate system
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On exit, a pointer to the created coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: NEW_COORDINATE_SYSTEM
    TYPE(COORDINATE_SYSTEM_PTR_TYPE), POINTER :: NEW_COORDINATE_SYSTEMS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_COORDINATE_SYSTEM)
    NULLIFY(NEW_COORDINATE_SYSTEMS)

    CALL ENTERS("COORDINATE_SYSTEM_CREATE_START",ERR,ERROR,*998)

    NULLIFY(NEW_COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,NEW_COORDINATE_SYSTEM,ERR,ERROR,*999)
    IF(ASSOCIATED(NEW_COORDINATE_SYSTEM)) THEN
      LOCAL_ERROR="Coordinate system number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
        & " has already been created."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
    ELSE
      IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
        CALL FLAG_ERROR("Coordinate system is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(NEW_COORDINATE_SYSTEM)
        ALLOCATE(NEW_COORDINATE_SYSTEM,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordinate system.",ERR,ERROR,*999)
      
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
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordinate systems.",ERR,ERROR,*999)
        DO coord_system_idx=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
          NEW_COORDINATE_SYSTEMS(coord_system_idx)%PTR=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR
        ENDDO !coord_system_idx
        NEW_COORDINATE_SYSTEMS(COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS+1)%PTR=>NEW_COORDINATE_SYSTEM
        DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
        COORDINATE_SYSTEMS%COORDINATE_SYSTEMS=>NEW_COORDINATE_SYSTEMS
        COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS+1
        
        COORDINATE_SYSTEM=>NEW_COORDINATE_SYSTEM
      ENDIF
    ENDIF
        
    CALL EXITS("COORDINATE_SYSTEM_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_COORDINATE_SYSTEM)) DEALLOCATE(NEW_COORDINATE_SYSTEM)
    IF(ASSOCIATED(NEW_COORDINATE_SYSTEMS)) DEALLOCATE(NEW_COORDINATE_SYSTEMS)
    NULLIFY(COORDINATE_SYSTEM)
998 CALL ERRORS("COORDINATE_SYSTEM_CREATE_START",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_CREATE_START")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a new coordinate system.
  SUBROUTINE COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to finish
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_idx

    CALL ENTERS("COORDINATE_SYSTEM_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
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

  !>Destroys a coordinate system.
  SUBROUTINE COORDINATE_SYSTEM_DESTROY(COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_no,new_coord_system_no
    LOGICAL :: FOUND
    TYPE(COORDINATE_SYSTEM_PTR_TYPE), POINTER :: NEW_COORDINATE_SYSTEMS(:)

    CALL ENTERS("COORDINATE_SYSTEM_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%USER_NUMBER==0) THEN
        CALL FLAG_ERROR("Cannot destroy the world coordinate system.",ERR,ERROR,*999)
      ELSE
        FOUND=.FALSE.
        new_coord_system_no=0
        ALLOCATE(NEW_COORDINATE_SYSTEMS(COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS-1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coordianate systems.",ERR,ERROR,*999)
        DO coord_system_no=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
          IF(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_no)%PTR%USER_NUMBER==COORDINATE_SYSTEM%USER_NUMBER) THEN
            FOUND=.TRUE.
          ELSE
            new_coord_system_no=new_coord_system_no+1
            NEW_COORDINATE_SYSTEMS(new_coord_system_no)%PTR=>COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_no)%PTR
          ENDIF
        ENDDO !coord_system_no
        IF(FOUND) THEN
          CALL COORDINATE_SYSTEM_FINALISE(COORDINATE_SYSTEM,ERR,ERROR,*999)
          DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
          COORDINATE_SYSTEMS%COORDINATE_SYSTEMS=>NEW_COORDINATE_SYSTEMS
          COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS-1
        ELSE
          DEALLOCATE(NEW_COORDINATE_SYSTEMS)
          CALL FLAG_ERROR("Coordinate system number to destroy does not exist.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("COORDINATE_SYSTEM_DESTROY")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEM_DESTROY",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEM_DESTROY")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DESTROY

  !
  !================================================================================================================================
  !
  
  !>Calculates DX(:)/DZ(I) at X, where Z(I) are rectangular cartesian and X(:) are curvilinear coordinates defined by COORDINATE_SYSTEM for double precision coordinates.
  FUNCTION DXZ_DP(COORDINATE_SYSTEM,I,X,ERR,ERROR)
  
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
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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

  !>Calculates the norm of a derivative in a coordinate system identified by COORDINATE_SYSTEM at the given interpolated
  !>point and returns the value in NORM for single precision coordinates. PART_DERIV_INDEX is used to select the
  !>appropriate partial derivative (i.e., wrt S1, S2 or S3) to normalise.
  SUBROUTINE COORDINATE_DERIVATIVE_NORM(COORDINATE_SYSTEM,PART_DERIV_INDEX,INTERPOLATED_POINT,DERIV_NORM,ERR,ERROR,*)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to calculate the derivative norm for
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_INDEX !<The partial derivative index to select the direction to normalise
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<A pointer to the interpolated point 
    REAL(DP), INTENT(OUT) :: DERIV_NORM !<On exit, the derivative norm of the coordinate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
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
              CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
                IF(NUMBER_OF_COMPONENTS==2) THEN
                  DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                    & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2
                ELSE IF(NUMBER_OF_COMPONENTS==3) THEN
                  DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                    & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2+INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX)**2
                ELSE
                  LOCAL_ERROR="The number of components for the interpolated point of "// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                    & " is invalid for a cylindrical polar coordinate system."
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
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              DERIV_NORM=SQRT(DERIV_NORM)
            CASE DEFAULT
              LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PART_DERIV_INDEX,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PART_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The interpolated point has a maximum number of partial derivatives of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATED_POINT%MAX_PARTIAL_DERIVATIVE_INDEX,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The point has not been interpolated to include first derivative values.",ERR,ERROR,*999)
        ENDIF          
      ELSE
        CALL FLAG_ERROR("Interpolated point is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
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

  !>Adjusts the interpolation for non-rectangular cartesian coordinate systems.
  SUBROUTINE COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,PARTIAL_DERIVATIVE_INDEX,VALUE,ERR,ERROR,*)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to adjust
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index to adjust
    REAL(DP), INTENT(INOUT) :: VALUE !<On entry, the coordinate value to adjust. On exit, the adjusted value.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: COSHX,CSS,D,DES,FOCUS,R,SS,SINHX,THETA
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("COORDINATE_INTERPOLATION_ADJUST",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        !Do nothing
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
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
            & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a cylindrical coordinate system."
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
            & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a cylindrical/spherical coordinate system."
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
            & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a prolate spheroidal coordinate system."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
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

  !>Adjusts the interpolation parameters for non-rectangular cartesian coordinate systems.
  SUBROUTINE COORDINATE_INTERPOLATION_PARAMETERS_ADJUST(COORDINATE_SYSTEM,INTERPOLATION_PARAMETERS,ERR,ERROR,*)
  
    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to adjust
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<A pointer to the interpolation parameters to adjust
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST",ERR,ERROR,*999)

!!TODO: Tidy up element parameters for non-rc coordinate systems. See bottom of XPXE and ZPZE.
    
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          !Do nothing
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a cylindrical coordinate system."
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
              & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a spherical coordinate system."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%RADIAL_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & RADIAL_INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for a prolate spheroidal coordinate system."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Interpolation parameters is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
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

!>Calculates transformation between spatial CS and rotated refernce orthogonal material CS.
  SUBROUTINE COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DXDNU, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT
    REAL(DP),INTENT(OUT) :: DXDNU(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS

    CALL ENTERS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE",ERR,ERROR,*999)

    NUMBER_OF_DIMENSIONS = GEOMETRIC_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS

    IF (NUMBER_OF_DIMENSIONS == 3) THEN
       CALL COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_3D(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DXDNU, &
          & ERR,ERROR,*999)
       ELSE
          CALL COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DXDNU, &
          & ERR,ERROR,*999)
       ENDIF

    CALL EXITS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE")
    RETURN
999 CALL ERRORS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE",ERR,ERROR)
    CALL EXITS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE")
    RETURN 1
  END SUBROUTINE

  !
  !================================================================================================================================
  !

!>Calculates transformation between spatial CS and rotated refernce orthogonal material CS in 2D space
  SUBROUTINE COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DXDNU, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT
    REAL(DP),INTENT(OUT) :: DXDNU(2,2)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx
    REAL(DP) :: ANGLE,DXDNUR(2,2),DXDXI(2,2),R(2,2),MAGNITUDE

    CALL ENTERS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D",ERR,ERROR,*999)

    derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1) !2,4,7
    DXDXI(1,1)=GEOMETRIC_INTERPOLATED_POINT%VALUES(1,derivative_idx) !dx1/dxi1
    DXDXI(2,1)=GEOMETRIC_INTERPOLATED_POINT%VALUES(2,derivative_idx) !dx2/dxi1

    !First calculate reference material CS
    !reference material direction 1.
    DXDNUR(:,1) = [ DXDXI(1,1),DXDXI(2,1) ]

    ! Compute (normalised) vector orthogonal to material direction 1 to form material direction 2
    DXDNUR(:,2) = [ -1*DXDNUR(2,1),DXDNUR(1,1) ]

    MAGNITUDE = L2NORM(DXDNUR(:,1))
    DXDNUR(1,1) = DXDNUR(1,1)/MAGNITUDE
    DXDNUR(2,1) = DXDNUR(2,1)/MAGNITUDE
    MAGNITUDE = L2NORM(DXDNUR(:,2))
    DXDNUR(1,2) = DXDNUR(1,2)/MAGNITUDE
    DXDNUR(2,2) = DXDNUR(2,2)/MAGNITUDE

    ANGLE = FIBRE_INTERPOLATED_POINT%VALUES(1,1)

    !Rotate by multiply with rotation matrix
    R(:,1) = [ COS(ANGLE),-SIN(ANGLE) ]
    R(:,2) = [ SIN(ANGLE),COS(ANGLE) ]

    CALL MATRIX_PRODUCT(R,DXDNUR,DXDNU,ERR,ERROR,*999)

    MAGNITUDE = L2NORM(DXDNU(:,1))
    DXDNU(1,1) = DXDNU(1,1)/MAGNITUDE
    DXDNU(2,1) = DXDNU(2,1)/MAGNITUDE

    MAGNITUDE = L2NORM(DXDNU(:,2))
    DXDNU(1,2) = DXDNU(1,2)/MAGNITUDE
    DXDNU(2,2) = DXDNU(2,2)/MAGNITUDE

    CALL EXITS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D")
    RETURN
999 CALL ERRORS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D",ERR,ERROR)
    CALL EXITS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D")
    RETURN 1
  END SUBROUTINE COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_2D

  !
  !================================================================================================================================
  !

!>Calculates transformation between spatial CS and rotated refernce orthogonal material CS in 3D space
  SUBROUTINE COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_3D(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DXDNU, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT
    REAL(DP),INTENT(OUT) :: DXDNU(3,3)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx,fibre_idx,geometric_idx,idx1,idx2,xi_idx
    INTEGER(INTG) :: NUMBER_OF_GEOMETRIC_COMPONENTS,NUMBER_OF_FIBRE_COMPONENTS,NUMBER_OF_XI_COORDS 
    INTEGER(INTG) :: vector(3) = (/1,2,3/)
    REAL(DP) :: ANGLE(3),DXDNU1(3,3),DXDNU2(3,3),DXDNU3(3,3),DXDNUR(3,3),DXDXI(3,3),F(3),G(3),H(3),Ra(3,3),Rb(3,3)
    
    CALL ENTERS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE",ERR,ERROR,*999)
    
    !initialse arrays
    DO idx1=1,3
      ANGLE(idx1)=0.0_DP
      DO idx2=1,3
        DXDNU(idx1,idx2)=0.0_DP
        DXDXI(idx1,idx2)=0.0_DP
      ENDDO
    ENDDO

    NUMBER_OF_GEOMETRIC_COMPONENTS=GEOMETRIC_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS 
    NUMBER_OF_FIBRE_COMPONENTS=FIBRE_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS    
    NUMBER_OF_XI_COORDS=GEOMETRIC_INTERPOLATED_POINT%interpolation_parameters%bases(1)%ptr%number_of_xi
    
    DO geometric_idx=1,NUMBER_OF_GEOMETRIC_COMPONENTS 
      DO xi_idx=1,NUMBER_OF_XI_COORDS
        derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7      
        DXDXI(geometric_idx,xi_idx)=GEOMETRIC_INTERPOLATED_POINT%VALUES(geometric_idx,derivative_idx) !dx/dxi
      ENDDO
    ENDDO

    DO fibre_idx=1,NUMBER_OF_FIBRE_COMPONENTS    
      ANGLE(fibre_idx)=FIBRE_INTERPOLATED_POINT%VALUES(fibre_idx,1) !fibre, imbrication and sheet. All in radians
    ENDDO

    !First calculate reference material CS
    DO geometric_idx=1,NUMBER_OF_GEOMETRIC_COMPONENTS
      F(geometric_idx)=DXDXI(geometric_idx,1) !reference material direction 1.
    ENDDO
    CALL CROSS_PRODUCT(DXDXI(vector,1),DXDXI(vector,2),H,ERR,ERROR,*999) !reference material direction 3.    
    CALL CROSS_PRODUCT(H,F,G,ERR,ERROR,*999) !reference material direction 2.      

    DO geometric_idx=1,NUMBER_OF_GEOMETRIC_COMPONENTS
      DXDNUR(geometric_idx,1)=F(geometric_idx)/L2NORM(F) 
      DXDNUR(geometric_idx,2)=G(geometric_idx)/L2NORM(G)      
      DXDNUR(geometric_idx,3)=H(geometric_idx)/L2NORM(H)        
    ENDDO

    !FIBRE ANGLE(alpha) - ANGLE(1) 
    !In order to rotate reference material CS by alpha(fibre angle) in anti-clockwise  
    !direction about its axis 3, following steps are performed.
    !(a) first align reference material direction 3 with Z(spatial) axis by rotating the ref matrial CS. 
    !(b) then rotate the aligned material CS by alpha about Z axis in anti-clockwise direction
    !(c) apply the inverse of step(a) to the CS in (b)
    !It can be shown that steps (a),(b) and (c) are equivalent to post-multiplying
    !rotation in (a) by rotation in (b). i.e. Ra*Rb  
       
    !The normalised reference material CS contains the transformation(rotaion) between 
    !the spatial CS -> reference material CS. i.e. Ra 
    DO idx1=1,3
      DO idx2=1,3
        Ra(idx1,idx2)=DXDNUR(idx1,idx2)  
      ENDDO
    ENDDO

    !Initialise rotation matrix Rb
    DO idx1=1,3
      DO idx2=1,3
        Rb(idx1,idx2)=0.0_DP 
        IF (idx1==idx2) THEN
          Rb(idx1,idx2)=1.0_DP  
        ENDIF
      ENDDO
    ENDDO    
    !Populate rotation matrix Rb about axis 3 (Z)
    Rb(1,1)=cos(ANGLE(1))
    Rb(1,2)=-sin(ANGLE(1))
    Rb(2,1)=sin(ANGLE(1))
    Rb(2,2)=cos(ANGLE(1))
    
    CALL MATRIX_PRODUCT(Ra,Rb,DXDNU3,ERR,ERROR,*999)  


    !IMBRICATION ANGLE (beta) - ANGLE(2)     
    !In order to rotate alpha-rotated material CS by beta(imbrication angle) in anti-clockwise  
    !direction about its new axis 2, following steps are performed.
    !(a) first align new material direction 2 with Y(spatial) axis by rotating the new matrial CS. 
    !(b) then rotate the aligned CS by beta about Y axis in anti-clockwise direction
    !(c) apply the inverse of step(a) to the CS in (b)
    !As mentioned above, (a),(b) and (c) are equivalent to post-multiplying
    !rotation in (a) by rotation in (b). i.e. Ra*Rb  
        
    !DXNU3 contains the transformation(rotaion) between 
    !the spatial CS -> alpha-rotated reference material CS. i.e. Ra 
    DO idx1=1,3
      DO idx2=1,3
        Ra(idx1,idx2)=DXDNU3(idx1,idx2)  
      ENDDO
    ENDDO   
    !Initialise rotation matrix Rb
    DO idx1=1,3
      DO idx2=1,3
        Rb(idx1,idx2)=0.0_DP 
        IF (idx1==idx2) THEN
          Rb(idx1,idx2)=1.0_DP  
        ENDIF
      ENDDO
    ENDDO    
    !Populate rotation matrix Rb about axis 2 (Y). Note the sign change
    Rb(1,1)=cos(ANGLE(2))
    Rb(1,3)=sin(ANGLE(2))
    Rb(3,1)=-sin(ANGLE(2))
    Rb(3,3)=cos(ANGLE(2))

    CALL MATRIX_PRODUCT(Ra,Rb,DXDNU2,ERR,ERROR,*999)  


    !SHEET ANGLE (gamma) - ANGLE(3)    
    !In order to rotate alpha-beta-rotated material CS by gama(sheet angle) in anti-clockwise  
    !direction about its new axis 1, following steps are performed.
    !(a) first align new material direction 1 with X(spatial) axis by rotating the new matrial CS. 
    !(b) then rotate the aligned CS by gama about X axis in anti-clockwise direction
    !(c) apply the inverse of step(a) to the CS in (b)
    !Again steps (a),(b) and (c) are equivalent to post-multiplying
    !rotation in (a) by rotation in (b). i.e. Ra*Rb  
        
    !DXNU2 contains the transformation(rotaion) between 
    !the spatial CS -> alpha-beta-rotated reference material CS. i.e. Ra 
    DO idx1=1,3
      DO idx2=1,3
        Ra(idx1,idx2)=DXDNU2(idx1,idx2)  
      ENDDO
    ENDDO   
    !Initialise rotation matrix Rb
    DO idx1=1,3
      DO idx2=1,3
        Rb(idx1,idx2)=0.0_DP 
        IF (idx1==idx2) THEN
          Rb(idx1,idx2)=1.0_DP  
        ENDIF
      ENDDO
    ENDDO    
    !Populate rotation matrix Rb about axis 1 (X). 
    Rb(2,2)=cos(ANGLE(3))
    Rb(2,3)=-sin(ANGLE(3))
    Rb(3,2)=sin(ANGLE(3))
    Rb(3,3)=cos(ANGLE(3))

    CALL MATRIX_PRODUCT(Ra,Rb,DXDNU1,ERR,ERROR,*999)  

    DO idx1=1,3
      DO idx2=1,3
        DXDNU(idx1,idx2)=DXDNU1(idx1,idx2)  
      ENDDO
    ENDDO   

    CALL EXITS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_3D")
    RETURN
999 CALL ERRORS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_3D",ERR,ERROR)
    CALL EXITS("COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_3D")
    RETURN 1
  END SUBROUTINE COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE_3D

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the coordinate system identified by USER_NUMBER. If a coordinate system with that number is not
  !>found then COORDINATE_SYSTEM is set to NULL.
  SUBROUTINE COORDINATE_SYSTEM_USER_NUMBER_FIND(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the coordinate system to find.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On exit, a pointer to the coordinate system with the specified user number if it exists. If no coordinate system has the specified user number the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    
    CALL ENTERS("COORDINATE_SYSTEM_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      CALL FLAG_ERROR("Coordinate_system is already associated.",ERR,ERROR,*999)
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

  !>Finalises the coordinate systems and destroys all coordinate systems.
  SUBROUTINE COORDINATE_SYSTEMS_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    
    CALL ENTERS("COORDINATE_SYSTEMS_FINALISE",ERR,ERROR,*999)

    DO coord_system_idx=1,COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS
      CALL COORDINATE_SYSTEM_FINALISE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(coord_system_idx)%PTR,ERR,ERROR,*999)
    ENDDO !coord_system_idx
    DEALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS)
    COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=0
    
    CALL EXITS("COORDINATE_SYSTEMS_FINALISE")
    RETURN
999 CALL ERRORS("COORDINATE_SYSTEMS_FINALISE",ERR,ERROR)
    CALL EXITS("COORDINATE_SYSTEMS_FINALISE")
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEMS_FINALISE

  !
  !================================================================================================================================
  !
   
  !>Initialises the coordinate systems and creates the world coordinate system.
  SUBROUTINE COORDINATE_SYSTEMS_INITIALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("COORDINATE_SYSTEMS_INITIALISE",ERR,ERROR,*999)

    !Allocate the coordinate systems
    ALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate coordinate systems.",ERR,ERROR,*999)
    !Create the default RC World cooordinate system
    ALLOCATE(COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR,STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate world coordinate system.",ERR,ERROR,*999)
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%USER_NUMBER=0
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%NUMBER_OF_DIMENSIONS=3
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%FOCUS=1.0_DP    
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%ORIGIN=(/0.0_DP,0.0_DP,0.0_DP/)
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%ORIENTATION=RESHAPE(&
      & (/1.0_DP,0.0_DP,0.0_DP, &
      &   0.0_DP,1.0_DP,0.0_DP, &
      &   0.0_DP,0.0_DP,1.0_DP/), &
      & (/3,3/))    
    COORDINATE_SYSTEMS%COORDINATE_SYSTEMS(1)%PTR%COORDINATE_SYSTEM_FINISHED=.TRUE.
    COORDINATE_SYSTEMS%NUMBER_OF_COORDINATE_SYSTEMS=1
   
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

