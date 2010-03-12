!> \file
!> $Id: region_routines.f90 690 2009-09-30 23:27:16Z chrispbradley $
!> \author David Nordsletten
!> \brief This module contains all region routines.
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

!> This module contains all region routines.
MODULE INTERFACE_ROUTINES

  USE BASE_ROUTINES
  USE COORDINATE_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE MESH_ROUTINES
  USE NODE_ROUTINES
!  USE REGION_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \brief Interface Type parameters
  INTEGER(INTG), PARAMETER :: INTERFACE_SURFACE_TYPE    =1
  INTEGER(INTG), PARAMETER :: INTERFACE_VOLUME_TYPE     =2
  INTEGER(INTG), PARAMETER :: INTERFACE_NODAL_TYPE      =1
  INTEGER(INTG), PARAMETER :: INTERFACE_INTEGRAL_TYPE   =2

  !Module types

  !Module variables

  TYPE(REGIONS_TYPE) :: REGIONS


  PUBLIC INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET, INTERFACE_EQUATIONS_SET_SETUP

  PUBLIC INTF_REGION_TYPE_SET ! <<>>

CONTAINS

  !
  !================================================================================================================================
  !


  !
  !================================================================================================================================
  !


  !>Finishes the creation of a region. \see OPENCMISS::CMISSRegionCreateFinish
  SUBROUTINE INTF_REGION_TYPE_SET(INTF_REGION,INTERFACE_TYPE,COUPLING_TYPE,ERR,ERROR,*)   ! <<>>

    !Argument variables
    TYPE(REGION_TYPE), POINTER ::       INTF_REGION    !<A pointer to the region to finish the creation of
    INTEGER(INTG), INTENT(IN) ::       INTERFACE_TYPE    !<The error code
    INTEGER(INTG), INTENT(IN) ::       COUPLING_TYPE    !<The error code
    INTEGER(INTG), INTENT(OUT) ::       ERR       !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR       !<The error string
    !Local Variables
    INTEGER(INTG) ::            I
     
    CALL ENTERS("INTF_REGION_TYPE_SET",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTF_REGION))             CALL FLAG_ERROR("Interface Region is not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTF_REGION%INTF))          CALL FLAG_ERROR("Interface region class is not associated.",ERR,ERROR,*999)
    !Setting the interface type (either nodal, or integral types)
    SELECT CASE(INTERFACE_TYPE)
    CASE(INTERFACE_NODAL_TYPE)
        INTF_REGION%INTF%INTERFACE_TYPE=INTERFACE_NODAL_TYPE
    CASE(INTERFACE_INTEGRAL_TYPE)
        INTF_REGION%INTF%INTERFACE_TYPE=INTERFACE_INTEGRAL_TYPE
    CASE DEFAULT
        CALL FLAG_ERROR("Improper interface type selected.",ERR,ERROR,*999)
    END SELECT
    !Setting the interface coupling type (either surface to surface, or volume)
    SELECT CASE(COUPLING_TYPE)
    CASE(INTERFACE_SURFACE_TYPE)
        INTF_REGION%INTF%COUPLING_TYPE=INTERFACE_SURFACE_TYPE
    CASE(INTERFACE_VOLUME_TYPE)
        INTF_REGION%INTF%COUPLING_TYPE=INTERFACE_VOLUME_TYPE
    CASE DEFAULT
        CALL FLAG_ERROR("Improper interface coupling type selected.",ERR,ERROR,*999)
    END SELECT

    CALL EXITS("INTF_REGION_TYPE_SET")
    RETURN
999 CALL ERRORS("INTF_REGION_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTF_REGION_TYPE_SET")
    RETURN 1
  END SUBROUTINE INTF_REGION_TYPE_SET  ! <<>>


  SUBROUTINE INTERFACE_MAPPING_CREATE_START(REGION,MESH,DOMAIN_MESH1,DOMAIN_MESH2,ERR,ERROR,*)
    
    !Argument variables
    TYPE(REGION_TYPE), POINTER ::       REGION !<A pointer to the region to create the mesh on
    TYPE(MESH_TYPE), POINTER ::       MESH !<A pointer to the interface mesh. Must be associated on entry.
    TYPE(MESH_TYPE), POINTER ::       DOMAIN_MESH1 !<A pointer to a domain mesh. Must be finished on entry.
    TYPE(MESH_TYPE), POINTER ::       DOMAIN_MESH2 !<A pointer to a domain mesh. Must be finished on entry.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::             I, J !<Dummy variables

    CALL ENTERS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR,*999)

    !Checking for potential errors on entry.
    IF(.NOT.ASSOCIATED(REGION))          CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(REGION%COORDINATE_SYSTEM))    CALL FLAG_ERROR("Coordinate system on region is not associated",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(REGION%INTF))         CALL FLAG_ERROR("Interface region component is not associated",ERR,ERROR,*999)
    IF(.NOT.REGION%INTERFACE_REGION)         CALL FLAG_ERROR("Region is not an interface region",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(MESH))            CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    IF(.NOT.MESH%MESH_FINISHED)            CALL FLAG_ERROR("Interface mesh is not finished.",ERR,ERROR,*999)
    IF(ASSOCIATED(MESH%INTF))            CALL FLAG_ERROR("Mesh interface map is already associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(DOMAIN_MESH1))          CALL FLAG_ERROR("First input domain mesh is not associated.",ERR,ERROR,*999)
    IF(.NOT.DOMAIN_MESH1%MESH_FINISHED)         CALL FLAG_ERROR("First input domain mesh is not finished.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(DOMAIN_MESH2))          CALL FLAG_ERROR("Second input domain mesh is not associated.",ERR,ERROR,*999)
    IF(.NOT.DOMAIN_MESH2%MESH_FINISHED)         CALL FLAG_ERROR("Second input domain mesh is not finished.",ERR,ERROR,*999)

    !Create the interface component to the mesh.
    ALLOCATE(MESH%INTF,STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface map.",ERR,ERROR,*999)
    !Assign pointer to the interface mesh.
    MESH%INTF%MESH=>MESH
    MESH%INTF%MAPPING_FINISHED=.FALSE.
    !Assign the number of coupled regions to the interface mesh and allocate the domain mapping.
    MESH%INTF%N_COUPLED_REGIONS=REGION%INTF%N_COUPLED_REGIONS
    ALLOCATE(MESH%INTF%DOMAIN_MAP(MESH%INTF%N_COUPLED_REGIONS))
    DO I = 1,MESH%INTF%N_COUPLED_REGIONS
       !Nullify and allocate the domain mapping pointer
       NULLIFY(MESH%INTF%DOMAIN_MAP(I)%PTR)
       ALLOCATE(MESH%INTF%DOMAIN_MAP(I)%PTR,STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface domain map.",ERR,ERROR,*999)
       !Nullify the domain mapping domain mesh which refers to the mesh which is to be mapped to the interface.
       MESH%INTF%DOMAIN_MAP(I)%PTR%GEOMETRIC_COMPONENT=0
       MESH%INTF%DOMAIN_MAP(I)%PTR%FIELD_COMPONENT=0
       !Declaring which mesh of region I is paired to the interface.
       NULLIFY(MESH%INTF%DOMAIN_MAP(I)%PTR%DOMAIN_MESH)
       IF(REGION%INTF%COUPLED_REGIONS(I)%PTR%USER_NUMBER==DOMAIN_MESH1%REGION%USER_NUMBER) THEN
         MESH%INTF%DOMAIN_MAP(I)%PTR%DOMAIN_MESH=>DOMAIN_MESH1
       ELSE IF (REGION%INTF%COUPLED_REGIONS(I)%PTR%USER_NUMBER==DOMAIN_MESH2%REGION%USER_NUMBER) THEN
         MESH%INTF%DOMAIN_MAP(I)%PTR%DOMAIN_MESH=>DOMAIN_MESH2
       ELSE
         CALL FLAG_ERROR("Neither input domain mesh is referenced to a region paired to the interface region.",ERR,ERROR,*999)
       END IF
       !Create the interface to domain map which defines for each interface element what domain elements it is connected to.
       IF(MESH%NUMBER_OF_ELEMENTS<=0) CALL FLAG_ERROR("Number of mesh elements must be > 0.",ERR,ERROR,*999)
       ALLOCATE(MESH%INTF%DOMAIN_MAP(I)%PTR%INTERFACE_TO_DOMAIN(MESH%NUMBER_OF_ELEMENTS))
       DO J = 1,MESH%NUMBER_OF_ELEMENTS
          !Allocate the interface to domain pointer and set the number of elements (domain elements) which are mapped to interface elements to zero.
          NULLIFY(MESH%INTF%DOMAIN_MAP(I)%PTR%INTERFACE_TO_DOMAIN(J)%PTR)
          ALLOCATE(MESH%INTF%DOMAIN_MAP(I)%PTR%INTERFACE_TO_DOMAIN(J)%PTR,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface domain, interface-to-domain map.",ERR,ERROR,*999)
          MESH%INTF%DOMAIN_MAP(I)%PTR%INTERFACE_TO_DOMAIN(J)%PTR%NUMBER_OF_ELEM=0
       END DO
    END DO

    CALL EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE INTERFACE_MAPPING_ELEMENT_ADD(MESH,MESH_ELEMENT,DOMAIN_MESH,DOMAIN_MESH_ELEMENTS,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_TYPE), POINTER ::       MESH !<A pointer to the interface mesh. Must be associated on entry.
    INTEGER(INTG) ::            MESH_ELEMENT !<Global element number for the interface element to be mapped to domain element list
    TYPE(MESH_TYPE), POINTER ::         DOMAIN_MESH !<Domain mesh.
    INTEGER(INTG), INTENT(IN) ::      DOMAIN_MESH_ELEMENTS(:) !<Elements in the domain mesh which are to be mapped to the interface element.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::             I !<Dummy variables
    INTEGER(INTG) ::             REGID !< Region ID variable

    CALL ENTERS("INTERFACE_MAPPING_ELEMENT_ADD",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(MESH))            CALL FLAG_ERROR("Interface mesh is not associated.",ERR,ERROR,*999)
    IF(.NOT.MESH%MESH_FINISHED)            CALL FLAG_ERROR("Interface mesh is not finished.",ERR,ERROR,*999)
    IF(MESH%INTF%MAPPING_FINISHED)         CALL FLAG_ERROR("Interface mapping is already finished.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(DOMAIN_MESH))          CALL FLAG_ERROR("Input domain mesh is not associated.",ERR,ERROR,*999)
    IF(.NOT.DOMAIN_MESH%MESH_FINISHED)         CALL FLAG_ERROR("Input domain mesh is not finished.",ERR,ERROR,*999)
    IF((MESH_ELEMENT<1).OR.(MESH_ELEMENT>MESH%NUMBER_OF_ELEMENTS)) THEN
       CALL FLAG_ERROR("Input domain mesh element does not exist.",ERR,ERROR,*999)
    END IF
    DO I=1,SIZE(DOMAIN_MESH_ELEMENTS,1)
       IF((DOMAIN_MESH_ELEMENTS(I)<1).OR.(DOMAIN_MESH_ELEMENTS(I)>DOMAIN_MESH%NUMBER_OF_ELEMENTS)) THEN
          CALL FLAG_ERROR("Input domain mesh and element list are incompatible.",ERR,ERROR,*999)
       END IF
    ENDDO
    REGID=0
    DO I=1,MESH%INTF%N_COUPLED_REGIONS
       IF(MESH%INTF%DOMAIN_MAP(I)%PTR%DOMAIN_MESH%REGION%USER_NUMBER==DOMAIN_MESH%REGION%USER_NUMBER) REGID=I
    ENDDO
    IF(REGID==0) CALL FLAG_ERROR("Input domain mesh is not associated with this interface.",ERR,ERROR,*999)
    IF(MESH%INTF%DOMAIN_MAP(REGID)%PTR%INTERFACE_TO_DOMAIN(MESH_ELEMENT)%PTR%NUMBER_OF_ELEM>0) THEN
      CALL FLAG_ERROR("Interface domain element is already finished.",ERR,ERROR,*999)
    END IF
    MESH%INTF%DOMAIN_MAP(REGID)%PTR%INTERFACE_TO_DOMAIN(MESH_ELEMENT)%PTR%NUMBER_OF_ELEM=SIZE(DOMAIN_MESH_ELEMENTS,1)
    ALLOCATE(MESH%INTF%DOMAIN_MAP(REGID)%PTR%INTERFACE_TO_DOMAIN(MESH_ELEMENT)%PTR%MAP(SIZE(DOMAIN_MESH_ELEMENTS,1)))
    MESH%INTF%DOMAIN_MAP(REGID)%PTR%INTERFACE_TO_DOMAIN(MESH_ELEMENT)%PTR%MAP=DOMAIN_MESH_ELEMENTS

    CALL EXITS("INTERFACE_MAPPING_ELEMENT_ADD")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_ELEMENT_ADD",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_ELEMENT_ADD")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH(MESH,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_TYPE), POINTER ::       MESH !<A pointer to the interface mesh. Must be associated on entry.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::             I, J !<Dummy variables

    CALL ENTERS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR,*999)

    !Checking for potential errors on entry.
    IF(.NOT.ASSOCIATED(MESH))            CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    IF(.NOT.MESH%MESH_FINISHED)            CALL FLAG_ERROR("Interface mesh is not finished.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(MESH%INTF))         CALL FLAG_ERROR("Mesh interface map is not associated.",ERR,ERROR,*999)
    IF(MESH%INTF%MAPPING_FINISHED)         CALL FLAG_ERROR("Mesh interface map is already finished.",ERR,ERROR,*999)
    MESH%INTF%MAPPING_FINISHED=.TRUE.

    CALL EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !


  SUBROUTINE INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_TYPE !<The equation type
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SUBTYPE !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_TYPE)
    CASE(EQUATIONS_SET_INTEGRAL_EQUATION_TYPE)
      CALL INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SUBTYPE,ERR,ERROR,*999)
    CASE(EQUATIONS_SET_NODAL_EQUATION_TYPE)
      CALL INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SUBTYPE,ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
        & " is not valid for a interface equations set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
       
    CALL EXITS("INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_SET_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !


  !>Sets/changes the equation subtype for a Laplace equation type of a classical field equations set class.
  SUBROUTINE INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    CASE(EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE)
      EQUATIONS_SET%CLASS=EQUATIONS_SET_INTERFACE_CLASS
      EQUATIONS_SET%TYPE=EQUATIONS_SET_INTEGRAL_EQUATION_TYPE
      EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE
    CASE(EQUATIONS_SET_NORMAL_KINEMATIC_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for an interface equations set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
       
    CALL EXITS("INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_INTEGRAL_EQUATIONS_SET_SUBTYPE_SET


  !
  !================================================================================================================================
  !


  !>Sets/changes the equation subtype for a Laplace equation type of a classical field equations set class.
  SUBROUTINE INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    CASE(EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE)
      EQUATIONS_SET%CLASS=EQUATIONS_SET_INTERFACE_CLASS
      EQUATIONS_SET%TYPE=EQUATIONS_SET_NODAL_EQUATION_TYPE
      EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE
    CASE(EQUATIONS_SET_NORMAL_KINEMATIC_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for an interface equations set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
       
    CALL EXITS("INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_NODAL_EQUATIONS_SET_SUBTYPE_SET


  !
  !================================================================================================================================
  !

  !
  SUBROUTINE INTERFACE_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_EQUATIONS_SET_SETUP",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_SET%TYPE)
    CASE(EQUATIONS_SET_INTEGRAL_EQUATION_TYPE)
      CALL INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
    CASE(EQUATIONS_SET_NODAL_EQUATION_TYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
        & " is not valid for a classical field equation set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
       
    CALL EXITS("INTERFACE_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Laplace equation type of a classical field equations set class.
  SUBROUTINE INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)

    SELECT CASE(EQUATIONS_SET%SUBTYPE)
    CASE(EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE)
      CALL INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
    CASE(EQUATIONS_SET_NORMAL_KINEMATIC_SUBTYPE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for a Laplace equation type of a classical field equation set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
       
    CALL EXITS("INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE INTERFACE_INTEGRAL_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the standard Laplace equation.
  SUBROUTINE INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS
    INTEGER(INTG):: DEPENDENT_FIELD_NUMBER_OF_VARIABLES,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG):: I
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    IF(EQUATIONS_SET%SUBTYPE /= EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE) THEN
       CALL FLAG_ERROR("Inappropriate equations set subtype.",ERR,ERROR,*999)
    END IF


    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, & 
          & ERR,ERROR,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION) ! DO Nothing
        CASE DEFAULT
           LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
           & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
           & " is invalid."
           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION) !Do nothing
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION) !Do nothing
        CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE) !Do nothing
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! DEPENDENT FIELD !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
         IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
            CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)  ! Initiate the build of the dependent field
            ! Initialize the field type
            CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999) 
            !define new created field to be dependent
            CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                                 & FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
            !start creation of a new field
            CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
            !look for decomposition rule already defined
            CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
            !apply decomposition rule found on new created field
            CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                                     & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
            !point new field to geometric field
            CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                                                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
            !set number of variables to 2 (1 for U and one for DELUDELN)
            DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
            CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                                      & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                                                 & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
            CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                            & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
            CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                            & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
            CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                            & FIELD_DP_TYPE,ERR,ERROR,*999)
            CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                            & FIELD_DP_TYPE,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE) THEN
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
            ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NORMAL_KINEMATIC_SUBTYPE) THEN
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
            END IF
            !calculate number of components with one component for each dimension
            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                                                       & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                                       & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            !set first i components to lambda or lambda . n fields


            STOP
!            DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
!               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I,<<>>, & 
!                                                     & ERR,ERROR,*999)
!               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
!                                                     & FIELD_DELUDELN_VARIABLE_TYPE,I,<<>>,ERR,ERROR,*999)  
               ! need to sort out how to flag the right field
!            END DO
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            !Specify fem solution method
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                           & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                           & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
              END DO
              CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solution method of " &
                         & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT  
         ELSE
            !Check the user specified field
            CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
            CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE, &
                                          & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
            CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                                     & ERR,ERROR,*999)
            CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                     & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
            CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
            CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                                     & ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            !calculate number of components with one component for each dimension and one for pressure
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_FULL_KINEMATIC_SUBTYPE) THEN
               DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
            ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NORMAL_KINEMATIC_SUBTYPE) THEN
               DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
            END IF
            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
               CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                                                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
               CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                                                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
            CASE DEFAULT
               LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                          &"*",ERR,ERROR))//" is invalid."
               CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
         END IF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
         IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) CALL FIELD_CREATE_FINISH( & 
            & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
      CASE DEFAULT
         LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for an interface"
         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Equation set setup equations type   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
       SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
       CASE(EQUATIONS_SET_SETUP_START_ACTION)
           IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
             CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
             CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
             CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
           ELSE
             CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
           ENDIF
       CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
             !Finish the equations creation
             CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
             CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
             !Create the equations mapping.
             CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
             CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,3,ERR,ERROR,*999)
             CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)   ! <<>>  I don't understand what's going on in this routine.
             CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
             CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
             CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
             SELECT CASE(EQUATIONS%SPARSITY_TYPE)
             CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                                                              & ERR,ERROR,*999)
             CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                                                              & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                                                                & ERR,ERROR,*999)
             CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
             END SELECT
             CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                     & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                     & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
       END SELECT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Equation set setup materials type   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
      CASE DEFAULT
      END SELECT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Equation set setup boundary conditions type   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE(EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE)
      SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
      CASE DEFAULT
      END SELECT  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Equation set setup analytic conditions type   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
      CASE DEFAULT
      END SELECT  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Equation set setup case default type   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE DEFAULT
       LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for standard interface set up."
       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
       
    CALL EXITS("INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Laplace equation type of an classical field equations set class.
  SUBROUTINE INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)

    SELECT CASE(SOLUTION_METHOD)
    CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
      EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
    CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    CALL EXITS("INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET
























END MODULE INTERFACE_ROUTINES