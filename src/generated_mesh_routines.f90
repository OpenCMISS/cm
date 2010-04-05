!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all generated mesh routines.
!> \todo Move generated regular mesh from mesh routines and generalise.
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
!> Contributor(s): Chris Bradley, Jack Lee
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

!> This module handles all generated mesh routines.
MODULE GENERATED_MESH_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshTypes GENERATED_MESH_ROUTINES::GeneratedMeshTypes 
  !> \brief Generated mesh types.
  !> \see GENERATED_MESH_ROUTINES,OPENCMISS_GeneratedMeshTypes
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_MESH_TYPE=1 !<A regular generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_POLAR_MESH_TYPE=2 !<A polar generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_FRACTAL_TREE_MESH_TYPE=3 !<A fractal tree generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_MESH_TYPE=4 !<A cylinder generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  !>@}

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces GENERATED_MESH_ROUTINES::GeneratedMeshCylinderSurfaces
  !> \brief Generated mesh cylinder type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_INNER_SURFACE=1  !<Inner surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_OUTER_SURFACE=2  !<Outer surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_TOP_SURFACE=3    !<Top surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_BOTTOM_SURFACE=4 !<Bottom surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  !>@}
  
  !Module types

  !Interfaces

  !>Starts the process of creating a generated mesh
  INTERFACE GENERATED_MESH_CREATE_START
    MODULE PROCEDURE GENERATED_MESH_CREATE_START_INTERFACE
    MODULE PROCEDURE GENERATED_MESH_CREATE_START_REGION
  END INTERFACE !GENERATED_MESH_CREATE_START

  !>Initialises the generated meshes for a region or interface.
  INTERFACE GENERATED_MESHES_INITIALISE
    MODULE PROCEDURE GENERATED_MESHES_INITIALISE_INTERFACE
    MODULE PROCEDURE GENERATED_MESHES_INITIALISE_REGION
  END INTERFACE !GENERATED_MESHES_INITIALISE

  !>Finds a generated mesh in a list of generated meshes in a region or interface.
  INTERFACE GENERATED_MESH_USER_NUMBER_FIND
    MODULE PROCEDURE GENERATED_MESH_USER_NUMBER_FIND_INTERFACE
    MODULE PROCEDURE GENERATED_MESH_USER_NUMBER_FIND_REGION
  END INTERFACE GENERATED_MESH_USER_NUMBER_FIND
  
  PUBLIC GENERATED_MESH_REGULAR_MESH_TYPE,GENERATED_MESH_POLAR_MESH_TYPE
  PUBLIC GENERATED_MESH_FRACTAL_TREE_MESH_TYPE,GENERATED_MESH_CYLINDER_MESH_TYPE

  public GENERATED_MESH_CYLINDER_INNER_SURFACE,GENERATED_MESH_CYLINDER_OUTER_SURFACE
  public GENERATED_MESH_CYLINDER_TOP_SURFACE,GENERATED_MESH_CYLINDER_BOTTOM_SURFACE
  
  PUBLIC GENERATED_MESHES_INITIALISE,GENERATED_MESHES_FINALISE

  PUBLIC GENERATED_MESH_CREATE_START,GENERATED_MESH_CREATE_FINISH

  PUBLIC GENERATED_MESH_DESTROY
  
  PUBLIC GENERATED_MESH_BASIS_SET,GENERATED_MESH_EXTENT_SET,GENERATED_MESH_NUMBER_OF_ELEMENTS_SET,GENERATED_MESH_ORIGIN_SET, &
    & GENERATED_MESH_TYPE_SET, GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE
    
  PUBLIC GENERATED_MESH_BASIS_GET,GENERATED_MESH_EXTENT_GET,GENERATED_MESH_NUMBER_OF_ELEMENTS_GET,GENERATED_MESH_ORIGIN_GET,&
    & GENERATED_MESH_TYPE_GET

  PUBLIC GENERATED_MESH_USER_NUMBER_FIND

  PUBLIC GENERATED_MESH_CYLINDER_SURFACES_GET
  
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Gets the basis of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBasisGet
  SUBROUTINE GENERATED_MESH_BASIS_GET(GENERATED_MESH,BASIS,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the basis of
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On return, the basis of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
            BASIS=>GENERATED_MESH%REGULAR_MESH%BASIS
          ELSE
            CALL FLAG_ERROR("Generated mesh regular mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)          
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            BASIS=>GENERATED_MESH%CYLINDER_MESH%BASIS
          ELSE
            CALL FLAG_ERROR("Generated mesh cylinder mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh generated type of "// &
            & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Generated mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_BASIS_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_BASIS_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_BASIS_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_BASIS_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBasisSet
  SUBROUTINE GENERATED_MESH_BASIS_SET(GENERATED_MESH,BASIS,ERR,ERROR,*)

     !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the basis of
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to generate the mesh with
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BASIS)) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
              GENERATED_MESH%REGULAR_MESH%BASIS=>BASIS 
            ELSE
              CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
              GENERATED_MESH%CYLINDER_MESH%BASIS=>BASIS 
            ELSE
              CALL FLAG_ERROR("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT          
        ELSE
          CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_BASIS_SET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_BASIS_SET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_BASIS_SET")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_BASIS_SET
  
  !
  !================================================================================================================================
  !
  
  !>Finishes the creation of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateFinish
  SUBROUTINE GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to finish the creation of
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<The mesh's user number
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the generated mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(MESH)) THEN
          CALL FLAG_ERROR("Mesh is already associated.",ERR,ERROR,*999)
        ELSE
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE) 
            CALL GENERATED_MESH_REGULAR_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Return the pointers
          MESH=>GENERATED_MESH%MESH
          MESH%GENERATED_MESH=>GENERATED_MESH
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_FINISH

 !
  !================================================================================================================================
  !

  !>Starts the creation of a generic generated mesh. 
  SUBROUTINE GENERATED_MESH_CREATE_START_GENERIC(GENERATED_MESHES,USER_NUMBER,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<A pointer to the generated meshes
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the generated mesh to create
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,generated_mesh_idx
    TYPE(GENERATED_MESH_TYPE), POINTER :: NEW_GENERATED_MESH
    TYPE(GENERATED_MESH_PTR_TYPE), POINTER :: NEW_GENERATED_MESHES(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    NULLIFY(NEW_GENERATED_MESH)
    NULLIFY(NEW_GENERATED_MESHES)

    CALL ENTERS("GENERATED_MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*997)
      ELSE
        !Initialise generated mesh
        CALL GENERATED_MESH_INITIALISE(NEW_GENERATED_MESH,ERR,ERROR,*999)
        !Set default generated mesh values
        NEW_GENERATED_MESH%USER_NUMBER=USER_NUMBER
        NEW_GENERATED_MESH%GLOBAL_NUMBER=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1
        NEW_GENERATED_MESH%GENERATED_MESHES=>GENERATED_MESHES        
        !Add new generated mesh into list of generated meshes
        ALLOCATE(NEW_GENERATED_MESHES(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new generated meshes.",ERR,ERROR,*999)
        DO generated_mesh_idx=1,GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES
          NEW_GENERATED_MESHES(generated_mesh_idx)%PTR=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
        ENDDO !generated_mesh_idx
        NEW_GENERATED_MESHES(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1)%PTR=>NEW_GENERATED_MESH
        IF(ASSOCIATED(GENERATED_MESHES%GENERATED_MESHES)) DEALLOCATE(GENERATED_MESHES%GENERATED_MESHES)
        GENERATED_MESHES%GENERATED_MESHES=>NEW_GENERATED_MESHES
        GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1
        !Return the pointer
        GENERATED_MESH=>NEW_GENERATED_MESH
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated meshes is not associated.",ERR,ERROR,*997)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_START_GENERIC")
    RETURN
999 CALL GENERATED_MESH_FINALISE(NEW_GENERATED_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_GENERATED_MESHES)) DEALLOCATE(NEW_GENERATED_MESHES)
    NULLIFY(GENERATED_MESH)
997 CALL ERRORS("GENERATED_MESH_CREATE_START_GENERIC",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_START_GENERIC")
    RETURN 1
    
  END SUBROUTINE GENERATED_MESH_CREATE_START_GENERIC
  
  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateFinish
  SUBROUTINE GENERATED_MESH_CREATE_START_INTERFACE(USER_NUMBER,INTERFACE,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the generated mesh to create
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the generated mesh on
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("GENERATED_MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(GENERATED_MESH)        
        CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,INTERFACE,GENERATED_MESH,ERR,ERROR,*999)
        IF(ASSOCIATED(GENERATED_MESH)) THEN
          LOCAL_ERROR="The specified user number of "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been used for a generated mesh."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(INTERFACE%GENERATED_MESHES)) THEN
            CALL GENERATED_MESH_CREATE_START_GENERIC(INTERFACE%GENERATED_MESHES,USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
            GENERATED_MESH%INTERFACE=>INTERFACE
          ELSE
            CALL FLAG_ERROR("Interface generated meshes is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_START_INTERFACE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_START_INTERFACE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_START_INTERFACE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateFinish
  SUBROUTINE GENERATED_MESH_CREATE_START_REGION(USER_NUMBER,REGION,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the generated mesh to create
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the generated mesh on
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("GENERATED_MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(GENERATED_MESH)
        CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,REGION,GENERATED_MESH,ERR,ERROR,*999)
        IF(ASSOCIATED(GENERATED_MESH)) THEN
          LOCAL_ERROR="The specified user number of "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been used for a generated mesh."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(REGION%GENERATED_MESHES)) THEN
            CALL GENERATED_MESH_CREATE_START_GENERIC(REGION%GENERATED_MESHES,USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
            GENERATED_MESH%REGION=>REGION
          ELSE
            CALL FLAG_ERROR("Region generated meshes is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_START_REGION")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_START_REGION",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_START_REGION")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh. \see OPENCMISS::CMISSGeneratedMeshCreateDestroy
  SUBROUTINE GENERATED_MESH_DESTROY(GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: generated_mesh_idx,generated_mesh_position
    TYPE(GENERATED_MESH_PTR_TYPE), POINTER :: NEW_GENERATED_MESHES(:)
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES

    CALL ENTERS("GENERATED_MESH_DESTROY",ERR,ERROR,*998)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      GENERATED_MESHES=>GENERATED_MESH%GENERATED_MESHES
      IF(ASSOCIATED(GENERATED_MESHES)) THEN
        IF(ASSOCIATED(GENERATED_MESHES%GENERATED_MESHES)) THEN
          generated_mesh_position=GENERATED_MESH%GLOBAL_NUMBER
          CALL GENERATED_MESH_FINALISE(GENERATED_MESH,ERR,ERROR,*999)        
          !Remove the generated mesh from the list of generated meshes
          IF(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES>1) THEN
            ALLOCATE(NEW_GENERATED_MESHES(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new generated meshes.",ERR,ERROR,*999)
            DO generated_mesh_idx=1,GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES
              IF(generated_mesh_idx<generated_mesh_position) THEN
                NEW_GENERATED_MESHES(generated_mesh_idx)%PTR=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
              ELSE IF(generated_mesh_idx>generated_mesh_position) THEN
                GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR%GLOBAL_NUMBER=GENERATED_MESHES% & 
                  & GENERATED_MESHES(generated_mesh_idx)%PTR%GLOBAL_NUMBER-1
                NEW_GENERATED_MESHES(generated_mesh_idx-1)%PTR=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
              ENDIF
            ENDDO !generated_mesh_idx
            DEALLOCATE(GENERATED_MESHES%GENERATED_MESHES)
            GENERATED_MESHES%GENERATED_MESHES=>NEW_GENERATED_MESHES
            GENERATED_MESHES%NUMBER_OF_GENERATED_MESHeS=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES-1
          ELSE
            DEALLOCATE(GENERATED_MESHES%GENERATED_MESHES)
            GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=0
          ENDIF
        ELSE
          CALL FLAG_ERROR("Generated meshes are not associated",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Generated mesh generated meshes is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*998)
    END IF
 
    CALL EXITS("GENERATED_MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_GENERATED_MESHES)) DEALLOCATE(NEW_GENERATED_MESHES)
998 CALL ERRORS("GENERATED_MESH_DESTROY",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_DESTROY")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Gets the extent of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshExtentGet
  SUBROUTINE GENERATED_MESH_EXTENT_GET(GENERATED_MESH,EXTENT,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: EXTENT(:) !<On return, maximum extent per axis, or inner & outer radii and length of cylinder
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_EXTENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(SIZE(EXTENT,1)>=SIZE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT,1)) THEN
          EXTENT=GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT
        ELSE
          LOCAL_ERROR="The size of EXTENT is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT,1),"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(SIZE(EXTENT,1)>=3) THEN
          EXTENT=GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT
        ELSE
          LOCAL_ERROR="The size of EXTENT is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))//" and it needs to be 3."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_EXTENT_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_EXTENT_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_EXTENT_GET")
    RETURN   
  END SUBROUTINE GENERATED_MESH_EXTENT_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the extent of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshExtentSet
  SUBROUTINE GENERATED_MESH_EXTENT_SET(GENERATED_MESH,EXTENT,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: EXTENT(:) !<The extent of the generated mesh (MAXIMUM for regular type, CYLINDER for cylinder type)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_EXTENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)) DEALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)
            ALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT(SIZE(EXTENT)),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate maximum extent.",ERR,ERROR,*999)
            GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT=EXTENT 
          ELSE
            CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN 
            ALLOCATE(GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT(SIZE(EXTENT)),STAT=ERR)
            GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT=EXTENT 
          ELSE
            CALL FLAG_ERROR("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_EXTENT_SET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_EXTENT_SET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_EXTENT_SET")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_EXTENT_SET

  !
  !================================================================================================================================
  !

  !>Finalises a generated mesh and dellocates all memory.
  SUBROUTINE GENERATED_MESH_FINALISE(GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%REGULAR_MESH,ERR,ERROR,*999)
      CALL GENERATED_MESH_CYLINDER_FINALISE(GENERATED_MESH%CYLINDER_MESH,ERR,ERROR,*999)
      NULLIFY(GENERATED_MESH%MESH%GENERATED_MESH)
      NULLIFY(GENERATED_MESH%MESH)
      DEALLOCATE(GENERATED_MESH)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_FINALISE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_FINALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_FINALISE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a generated mesh.
  SUBROUTINE GENERATED_MESH_INITIALISE(GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(GENERATED_MESH,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate generated mesh.",ERR,ERROR,*999)
      GENERATED_MESH%USER_NUMBER=0
      GENERATED_MESH%GLOBAL_NUMBER=0
      GENERATED_MESH%GENERATED_MESH_FINISHED=.FALSE.
      NULLIFY(GENERATED_MESH%REGION)
      NULLIFY(GENERATED_MESH%INTERFACE)
      GENERATED_MESH%GENERATED_TYPE=0
      NULLIFY(GENERATED_MESH%REGULAR_MESH)
      NULLIFY(GENERATED_MESH%CYLINDER_MESH)
      NULLIFY(GENERATED_MESH%MESH)
      !Default to a regular mesh.
      CALL GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_INITIALISE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_INITIALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_INITIALISE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the number of elements in a generated mesh.  \see OPENCMISS::CMISSGeneratedMeshNumberOfElementsGet
  SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_GET(GENERATED_MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ELEMENTS(:) !<On return, number of elements per axis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(SIZE(NUMBER_OF_ELEMENTS,1)>=SIZE(GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI,1)) THEN
          NUMBER_OF_ELEMENTS=GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
          LOCAL_ERROR="The size of NUMBER_OF_ELEMENTS is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(SIZE(NUMBER_OF_ELEMENTS,1)>=SIZE(GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI,1)) THEN
          NUMBER_OF_ELEMENTS=GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
          LOCAL_ERROR="The size of NUMBER_OF_ELEMENTS is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a generated mesh. \see OPENCMISS::CMISSGeneratedMeshNumberOfElementsSet
  SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,NUMBER_OF_ELEMENTS_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI(:) !<NUMBER_OF_ELEMENTS_XI(ni). The number of elements in the ni'th xi direction (or, r,theta,z for cylinder) to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
            ALLOCATE(GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(SIZE(NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
            GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI 
          ELSE
            CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            ALLOCATE(GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI(SIZE(NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
            GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI
          ELSE
            CALL FLAG_ERROR("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET
  
  !
  !================================================================================================================================
  !

  !>Get the origin of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshOriginGet
  SUBROUTINE GENERATED_MESH_ORIGIN_GET(GENERATED_MESH,ORIGIN,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: ORIGIN(:) !<On return, the origin coordinate for each axis  
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_ORIGIN_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(SIZE(ORIGIN,1)>=SIZE(GENERATED_MESH%REGULAR_MESH%ORIGIN,1)) THEN
          ORIGIN=GENERATED_MESH%REGULAR_MESH%ORIGIN 
        ELSE
          LOCAL_ERROR="The size of ORIGIN is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%REGULAR_MESH%ORIGIN,1),"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(SIZE(ORIGIN,1)>=SIZE(GENERATED_MESH%CYLINDER_MESH%ORIGIN,1)) THEN
          ORIGIN=GENERATED_MESH%CYLINDER_MESH%ORIGIN 
        ELSE
          LOCAL_ERROR="The size of ORIGIN is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))//" and it needs to be 3."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_ORIGIN_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ORIGIN_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ORIGIN_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ORIGIN_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshOriginSet
  SUBROUTINE GENERATED_MESH_ORIGIN_SET(GENERATED_MESH,ORIGIN,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: ORIGIN(:) !<ORIGIN(nj). The nj'th coordinate origin for the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_ORIGIN_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
            IF(.NOT.ALLOCATED(GENERATED_MESH%REGULAR_MESH%ORIGIN)) THEN
              ALLOCATE(GENERATED_MESH%REGULAR_MESH%ORIGIN(SIZE(ORIGIN)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
            ENDIF
            GENERATED_MESH%REGULAR_MESH%ORIGIN=ORIGIN 
          ELSE
            CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN 
            IF(SIZE(ORIGIN,1)==3) THEN
              IF(.NOT.ALLOCATED(GENERATED_MESH%CYLINDER_MESH%ORIGIN)) THEN
                ALLOCATE(GENERATED_MESH%CYLINDER_MESH%ORIGIN(SIZE(ORIGIN)),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Cylinder generated mesh is only supported for 3D.",ERR,ERROR,*999)
            ENDIF
            GENERATED_MESH%CYLINDER_MESH%ORIGIN=ORIGIN 
          ELSE
            CALL FLAG_ERROR("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
          END IF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_ORIGIN_SET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ORIGIN_SET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ORIGIN_SET")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_ORIGIN_SET
  
  !
  !================================================================================================================================
  ! 
  
  !>Start to create the regular generated mesh type
  SUBROUTINE GENERATED_MESH_REGULAR_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<The user number for the mesh to generate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: origin_idx
    INTEGER(INTG), ALLOCATABLE :: NUMBER_ELEMENTS_XI(:)
    INTEGER(INTG) :: ni,ne,ne1,ne2,ne3,NN,nn1,nn2,nn3,np,TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_ELEMENTS_XI(3), &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: REGION 
    
    CALL ENTERS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(GENERATED_MESH)) THEN
      REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH        
      IF(ASSOCIATED(REGULAR_MESH)) THEN
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          COORDINATE_SYSTEM=>REGION%COORDINATE_SYSTEM
        ELSE
          INTERFACE=>GENERATED_MESH%INTERFACE
          IF(ASSOCIATED(INTERFACE)) THEN
            IF(ASSOCIATED(INTERFACE%PARENT_REGION)) THEN
              COORDINATE_SYSTEM=>INTERFACE%PARENT_REGION%COORDINATE_SYSTEM
            ELSE
              CALL FLAG_ERROR("Generated mesh interface parent region is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Generated mesh does not have a region or interface associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          !TODO is regular type only for COORDINATE_RECTANGULAR_CARTESIAN_TYPE? 
          !If that, should we use IF rather than select?
          SELECT CASE(COORDINATE_SYSTEM%TYPE)
          CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
            !Determine the coordinate system and create the regular mesh for that system
            REGULAR_MESH%MESH_DIMENSION=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            NUMBER_OF_DIMENSIONS=REGULAR_MESH%MESH_DIMENSION
            IF(.NOT.ALLOCATED(REGULAR_MESH%ORIGIN)) THEN
              ALLOCATE(REGULAR_MESH%ORIGIN(NUMBER_OF_DIMENSIONS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
              DO origin_idx=1,NUMBER_OF_DIMENSIONS
                REGULAR_MESH%ORIGIN(origin_idx)=0.0_DP
              ENDDO !origin_idx
            ENDIF
            IF(SIZE(REGULAR_MESH%ORIGIN)==REGULAR_MESH%MESH_DIMENSION) THEN
              IF(SIZE(REGULAR_MESH%MAXIMUM_EXTENT)==REGULAR_MESH%MESH_DIMENSION) THEN
                IF(ASSOCIATED(REGULAR_MESH%BASIS)) THEN
                  BASIS=>REGULAR_MESH%BASIS
                  ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                  NUMBER_ELEMENTS_XI=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
                  SELECT CASE(BASIS%TYPE)
                  CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                    IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                      IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) CALL FLAG_ERROR("Must have 1 or more elements in all directions.", &
                        & ERR,ERROR,*999)
                      IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                        & CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
                      !Calculate sizes
                      TOTAL_NUMBER_OF_NODES=1
                      TOTAL_NUMBER_OF_ELEMENTS=1
                      TOTAL_NUMBER_OF_NODES_XI=1
                      TOTAL_NUMBER_ELEMENTS_XI=0
                      DO ni=1,BASIS%NUMBER_OF_XI
                        TOTAL_NUMBER_OF_NODES_XI(ni)=(BASIS%NUMBER_OF_NODES_XI(ni)-2)*NUMBER_ELEMENTS_XI(ni)+ &
                          & NUMBER_ELEMENTS_XI(ni)+1
                        TOTAL_NUMBER_ELEMENTS_XI(ni)=NUMBER_ELEMENTS_XI(ni)
                        TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES*TOTAL_NUMBER_OF_NODES_XI(ni)
                        TOTAL_NUMBER_OF_ELEMENTS=TOTAL_NUMBER_OF_ELEMENTS*TOTAL_NUMBER_ELEMENTS_XI(ni)
                      ENDDO !ni
                      !Create the default node set
                      !TODO we finish create after the nodes are initialised?
                      NULLIFY(NODES)
                      IF(ASSOCIATED(REGION)) THEN
                        CALL NODES_CREATE_START(REGION,TOTAL_NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
                      ELSE
                        CALL NODES_CREATE_START(INTERFACE,TOTAL_NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
                      ENDIF
                      !Finish the nodes creation
                      CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
                      !Create the mesh
                      IF(ASSOCIATED(REGION)) THEN
                        CALL MESH_CREATE_START(MESH_USER_NUMBER,REGION,SIZE(NUMBER_ELEMENTS_XI,1),GENERATED_MESH%MESH, &
                          & ERR,ERROR,*999)
                      ELSE
                        CALL MESH_CREATE_START(MESH_USER_NUMBER,INTERFACE,SIZE(NUMBER_ELEMENTS_XI,1),GENERATED_MESH%MESH, &
                          & ERR,ERROR,*999)
                      ENDIF
                      !Create the elements
                      CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
                      NULLIFY(MESH_ELEMENTS)
                      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,1,BASIS,MESH_ELEMENTS,ERR,ERROR,*999)
                      !Set the elements for the regular mesh
                      ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element nodes",ERR,ERROR,*999)                
                      !Step in the xi(3) direction
                      DO ne3=1,TOTAL_NUMBER_ELEMENTS_XI(3)+1
                        DO ne2=1,TOTAL_NUMBER_ELEMENTS_XI(2)+1
                          DO ne1=1,TOTAL_NUMBER_ELEMENTS_XI(1)+1
                            IF(BASIS%NUMBER_OF_XI<3.OR.ne3<=TOTAL_NUMBER_ELEMENTS_XI(3)) THEN
                              IF(BASIS%NUMBER_OF_XI<2.OR.ne2<=TOTAL_NUMBER_ELEMENTS_XI(2)) THEN
                                IF(ne1<=TOTAL_NUMBER_ELEMENTS_XI(1)) THEN
                                  ne=ne1
                                  np=1+(ne1-1)*(BASIS%NUMBER_OF_NODES_XI(1)-1)
                                  IF(BASIS%NUMBER_OF_XI>1) THEN
                                    ne=ne+(ne2-1)*TOTAL_NUMBER_ELEMENTS_XI(1)
                                    np=np+(ne2-1)*TOTAL_NUMBER_OF_NODES_XI(1)*(BASIS%NUMBER_OF_NODES_XI(2)-1)
                                    IF(BASIS%NUMBER_OF_XI>2) THEN
                                      ne=ne+(ne3-1)*TOTAL_NUMBER_ELEMENTS_XI(1)*TOTAL_NUMBER_ELEMENTS_XI(2)
                                      np=np+(ne3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)* &
                                        & (BASIS%NUMBER_OF_NODES_XI(3)-1)
                                    ENDIF
                                  ENDIF
                                  nn=0
                                  DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                    nn=nn+1
                                    ELEMENT_NODES(nn)=np+(nn1-1)                              
                                  ENDDO !nn1
                                  IF(BASIS%NUMBER_OF_XI>1) THEN
                                    DO nn2=2,BASIS%NUMBER_OF_NODES_XI(2)
                                      DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                        nn=nn+1
                                        ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)
                                      ENDDO !nn1
                                    ENDDO !nn2
                                    IF(BASIS%NUMBER_OF_XI>2) THEN
                                      DO nn3=2,BASIS%NUMBER_OF_NODES_XI(3)
                                        DO nn2=1,BASIS%NUMBER_OF_NODES_XI(2)
                                          DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                            nn=nn+1
                                            ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)+ &
                                              & (nn3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                          ENDDO !nn1
                                        ENDDO !nn2
                                      ENDDO !nn3
                                    ENDIF
                                  ENDIF
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS,ELEMENT_NODES, &
                                    & ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDIF
                          ENDDO !ne1
                        ENDDO !ne2
                      ENDDO !ne3
                      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,ERR,ERROR,*999)
                      !Finish the mesh
                      CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,ERR,ERROR,*999)                        
                    ELSE
                      CALL FLAG_ERROR("The number of xi directions of the given basis does not match the size of &
                        &the number of elements for the mesh.",ERR,ERROR,*999)
                    ENDIF
                  CASE(BASIS_SIMPLEX_TYPE)                  
                    CALL FLAG_ERROR("Regular meshes with simplex basis types is not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    CALL FLAG_ERROR("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                  &the maximum extent.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                &the origin.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            CALL FLAG_ERROR("Coordinate type is either invalid or not implemented.",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Coordiate System is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Regular mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    
    CALL EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    CALL ERRORS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_CREATE_FINISH
  
  !
  !================================================================================================================================
  ! 
  
  !>Start to create the regular generated mesh type
  SUBROUTINE GENERATED_MESH_CYLINDER_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), ALLOCATABLE :: NUMBER_ELEMENTS_XI(:)
    TYPE(REGION_TYPE), POINTER :: REGION 
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: ne1,ne2,ne3,nn1,nn2,nn3,from1,from2,from3,nn,ne
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    REAL(DP) :: DELTA(3),DELTAi(3)
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    
    CALL ENTERS("GENERATED_MESH_CYLINDER_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CYLINDER_MESH=>GENERATED_MESH%CYLINDER_MESH        
      IF(ASSOCIATED(CYLINDER_MESH)) THEN
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
            !TODO is regular type only for COORDINATE_RECTANGULAR_CARTESIAN_TYPE? 
            !If that, should we use IF rather than select?
            SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              !Determine the coordinate system and create the regular mesh for that system
              CYLINDER_MESH%MESH_DIMENSION=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              NUMBER_OF_DIMENSIONS=CYLINDER_MESH%MESH_DIMENSION
              IF(NUMBER_OF_DIMENSIONS==3) THEN ! hard-coded for 3D only
                IF(.NOT.ALLOCATED(CYLINDER_MESH%ORIGIN)) THEN
                  ALLOCATE(CYLINDER_MESH%ORIGIN(NUMBER_OF_DIMENSIONS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
                  CYLINDER_MESH%ORIGIN=0.0_DP
                ENDIF
                IF(SIZE(CYLINDER_MESH%ORIGIN)==CYLINDER_MESH%MESH_DIMENSION) THEN
                  IF(SIZE(CYLINDER_MESH%CYLINDER_EXTENT)==CYLINDER_MESH%MESH_DIMENSION) THEN
                    IF(ASSOCIATED(CYLINDER_MESH%BASIS)) THEN
                      BASIS=>CYLINDER_MESH%BASIS
                      ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                      NUMBER_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
                      SELECT CASE(BASIS%TYPE)
                      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                        IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                          IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) CALL FLAG_ERROR("Must have 1 or more elements in all directions.", &
                            & ERR,ERROR,*999)
                          IF(NUMBER_ELEMENTS_XI(2)<3) CALL FLAG_ERROR("Need >2 elements around the circumferential direction.", &
                            & ERR,ERROR,*999)
                          IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                            & CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
                          !Calculate nodes and element sizes
                          CALL GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,BASIS%NUMBER_OF_NODES_XI, &
                            & CYLINDER_MESH%CYLINDER_EXTENT, TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS, &
                            & NIDX,EIDX,DELTA,DELTAi,ERR,ERROR,*999)
                          !Create the default node set
                          !TODO we finish create after the nodes are initialised?
                          NULLIFY(NODES)
                          CALL NODES_CREATE_START(REGION,TOTAL_NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
                          !Finish the nodes creation
                          CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
                          !Create the mesh
                          CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%REGION,SIZE(NUMBER_ELEMENTS_XI,1), &
                            & GENERATED_MESH%MESH,ERR,ERROR,*999)
                          !Create the elements
                          CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
                          NULLIFY(MESH_ELEMENTS)
                          CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,1,BASIS,MESH_ELEMENTS,ERR,ERROR,*999)
                          !Set the elements for the regular mesh
                          ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element nodes.",ERR,ERROR,*999)                
                          ! calculate element topology (nodes per each element)
                          ! the idea is to translate given (r,theta,z) to NIDX equivalents, which include interior nodes
                          ne=0
                          DO ne3=1,NUMBER_ELEMENTS_XI(3)
                            from3=NINT(DELTA(3)*(ne3-1)/DELTAi(3)+1)
                            DO ne2=1,NUMBER_ELEMENTS_XI(2)
                              from2=NINT(DELTA(2)*(ne2-1)/DELTAi(2)+1)
                              DO ne1=1,NUMBER_ELEMENTS_XI(1)
                                from1=NINT(DELTA(1)*(ne1-1)/DELTAi(1)+1)
                                nn=0
                                ! number of nodes in an element is dependent on basis used
                                DO nn3=from3,from3+BASIS%NUMBER_OF_NODES_XI(3)-1
                                  DO nn2=from2,from2+BASIS%NUMBER_OF_NODES_XI(2)-1
                                    DO nn1=from1,from1+BASIS%NUMBER_OF_NODES_XI(1)-1
                                      nn=nn+1
                                      ! compensate for circumferential loop-around
                                      IF(nn2>SIZE(NIDX,2)) THEN
                                        ! DEBUG: little check here
                                        IF(nn2>SIZE(NIDX,2)+1) CALL FLAG_ERROR("NIDX needs debugging",ERR,ERROR,*999)
                                        ELEMENT_NODES(nn)=NIDX(nn1,1,nn3)
                                      ELSE
                                        ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                      ENDIF
                                    ENDDO ! nn1
                                  ENDDO ! nn2
                                ENDDO ! nn3
                                ne=ne+1
                                CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS,ELEMENT_NODES,ERR,ERROR,*999)
                              ENDDO ! ne1
                            ENDDO ! ne2
                          ENDDO ! ne3
                          CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,ERR,ERROR,*999)
                          !Finish the mesh
                          CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,ERR,ERROR,*999)                        
                        ELSE
                          CALL FLAG_ERROR("The number of xi directions of the given basis does not match the size of &
                            &the number of elements for the mesh.",ERR,ERROR,*999)
                        ENDIF  
                      CASE(BASIS_SIMPLEX_TYPE)                  
                        CALL FLAG_ERROR("Cylinder meshes with simplex basis types is not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        CALL FLAG_ERROR("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                      &the maximum extent.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                    &the origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Cylinder mesh requires a 3 dimensional coordinate system.",ERR,ERROR,*999) 
              ENDIF
            CASE DEFAULT
              CALL FLAG_ERROR("Coordinate type is either invalid or not implemented.",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Coordiate System is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Regular mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_CYLINDER_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    CALL ERRORS("GENERATED_MESH_CYLINDER_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_CREATE_FINISH

  !
  !================================================================================================================================
  ! 
  
  !>Finalise the cylinder mesh type
  SUBROUTINE GENERATED_MESH_CYLINDER_FINALISE(CYLINDER_MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("GENERATED_MESH_CYLINDER_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CYLINDER_MESH)) THEN
      IF(ALLOCATED(CYLINDER_MESH%ORIGIN)) DEALLOCATE(CYLINDER_MESH%ORIGIN)
      IF(ALLOCATED(CYLINDER_MESH%CYLINDER_EXTENT)) DEALLOCATE(CYLINDER_MESH%CYLINDER_EXTENT)
      IF(ALLOCATED(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)) DEALLOCATE(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)
      DEALLOCATE(CYLINDER_MESH)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_CYLINDER_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL ERRORS("GENERATED_MESH_CYLINDER_FINALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_FINALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_FINALISE

  !
  !================================================================================================================================
  ! 
  
  !>Initialise the regular generated mesh type
  SUBROUTINE GENERATED_MESH_CYLINDER_INITIALISE(GENERATED_MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("GENERATED_MESH_CYLINDER_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
        CALL FLAG_ERROR("Cylinder mesh is already associated for this generated mesh.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%CYLINDER_MESH,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate cylinder generated mesh.",ERR,ERROR,*999)
        GENERATED_MESH%CYLINDER_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_CYLINDER_MESH_TYPE
        NULLIFY(GENERATED_MESH%CYLINDER_MESH%BASIS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_CYLINDER_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_CYLINDER_FINALISE(GENERATED_MESH%CYLINDER_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GENERATED_MESH_CYLINDER_INITIALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_INITIALISE

  !
  !================================================================================================================================
  ! 
  
  !>Finalise the regular mesh type
  SUBROUTINE GENERATED_MESH_REGULAR_FINALISE(REGULAR_MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("GENERATED_MESH_REGULAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGULAR_MESH)) THEN
      IF(ALLOCATED(REGULAR_MESH%ORIGIN)) DEALLOCATE(REGULAR_MESH%ORIGIN)
      IF(ALLOCATED(REGULAR_MESH%MAXIMUM_EXTENT)) DEALLOCATE(REGULAR_MESH%MAXIMUM_EXTENT)
      IF(ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) DEALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)
      DEALLOCATE(REGULAR_MESH)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_REGULAR_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL ERRORS("GENERATED_MESH_REGULAR_FINALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_FINALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_FINALISE
  
  !
  !================================================================================================================================
  ! 
  
  !>Initialise the regular generated mesh type
  SUBROUTINE GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
     
    CALL ENTERS("GENERATED_MESH_REGULAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
        CALL FLAG_ERROR("Regular mesh is already associated for this generated mesh.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%REGULAR_MESH,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate regular generated mesh.",ERR,ERROR,*999)
        GENERATED_MESH%REGULAR_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_REGULAR_MESH_TYPE
        NULLIFY(GENERATED_MESH%REGULAR_MESH%BASIS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%REGULAR_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GENERATED_MESH_REGULAR_INITIALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the type of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshTypeGet
  SUBROUTINE GENERATED_MESH_TYPE_GET(GENERATED_MESH,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the type of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      TYPE=GENERATED_MESH%GENERATED_TYPE
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_TYPE_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_TYPE_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_TYPE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshTypeSet
  SUBROUTINE GENERATED_MESH_TYPE_SET(GENERATED_MESH,GENERATED_TYPE,ERR,ERROR,*)

     !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: GENERATED_TYPE !<The type of mesh to generate \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: OLD_GENERATED_TYPE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has already been finished.",ERR,ERROR,*999)
      ELSE
        OLD_GENERATED_TYPE=GENERATED_MESH%GENERATED_TYPE
        IF(OLD_GENERATED_TYPE/=GENERATED_TYPE) THEN
          !Initialise the new generated mesh type
          SELECT CASE(GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE) 
            CALL GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The specified generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Finalise the new generated mesh type          
          SELECT CASE(OLD_GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE) 
            CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%REGULAR_MESH,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_FINALISE(GENERATED_MESH%CYLINDER_MESH,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(OLD_GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_TYPE_SET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_TYPE_SET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_TYPE_SET")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_TYPE_SET
  
  !
  !================================================================================================================================
  !
  
  !>Finds and returns in generated mesh a pointer to that identified by USER_NUMBER in the given list of GENERATED_MESHES. If no generated mesh with that number exists GENERATED_MESH is left nullified.
  SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,GENERATED_MESHES,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<A pointer to the generated meshes to find the user number in
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: generated_mesh_idx

    CALL ENTERS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN      
        CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(GENERATED_MESH)
        generated_mesh_idx=1
        DO WHILE(generated_mesh_idx<=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES.AND..NOT.ASSOCIATED(GENERATED_MESH))
          IF(GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            GENERATED_MESH=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
            EXIT
          ELSE
            generated_mesh_idx=generated_mesh_idx+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated meshes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC")
    RETURN
999 CALL ERRORS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND_GENERIC")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND_GENERIC

  !
  !================================================================================================================================
  !
  
  !>Finds and returns in generated mesh a pointer to that identified by USER_NUMBER in the given INTERFACE. If no generated mesh with that number exists GENERATED MESH is left nullified.
  SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND_INTERFACE(USER_NUMBER,INTERFACE,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the generated mesh
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      CALL GENERATED_MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%GENERATED_MESHES,GENERATED_MESH,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND_INTERFACE

  !
  !================================================================================================================================
  !
  
  !>Finds and returns in generated mesh a pointer to that identified by USER_NUMBER in the given REGION. If no generated mesh with that number exists GENERATED MESH is left nullified.
  SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND_REGION(USER_NUMBER,REGION,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the generated mesh
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_USER_NUMBER_FIND_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL GENERATED_MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%GENERATED_MESHES,GENERATED_MESH,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND_REGION")
    RETURN
999 CALL ERRORS("GENERATED_MESH_USER_NUMBER_FIND_REGION",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND_REGION")
    RETURN 1
    
  END SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND_REGION

  !
  !================================================================================================================================
  !

  !>Finalises all generated meshes and deallocates all memory.
  SUBROUTINE GENERATED_MESHES_FINALISE(GENERATED_MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<A pointer to the generated meshes to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH

    CALL ENTERS("GENERATED_MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      DO WHILE(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES>0)
        GENERATED_MESH=>GENERATED_MESHES%GENERATED_MESHES(1)%PTR
        CALL GENERATED_MESH_DESTROY(GENERATED_MESH,ERR,ERROR,*999)
      ENDDO !generated_mesh_idx
      DEALLOCATE(GENERATED_MESHES)
    ENDIF
    
    CALL EXITS("GENERATED_MESHES_FINALISE")
    RETURN
999 CALL ERRORS("GENERATED_MESHES_FINALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESHES_FINALISE")
    RETURN 1
    
  END SUBROUTINE GENERATED_MESHES_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all generated meshes.
  SUBROUTINE GENERATED_MESHES_INITIALISE_GENERIC(GENERATED_MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESHES_TYPE), POINTER :: GENERATED_MESHES !<A pointer to the generated meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("GENERATED_MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(GENERATED_MESHES)) THEN
      CALL FLAG_ERROR("Generated meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(GENERATED_MESHES,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Generated meshes is not associated.",ERR,ERROR,*999)
      GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=0
      NULLIFY(GENERATED_MESHES%GENERATED_MESHES)
      NULLIFY(GENERATED_MESHES%REGION)
      NULLIFY(GENERATED_MESHES%INTERFACE)
    ENDIF
      
    CALL EXITS("GENERATED_MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL GENERATED_MESHES_FINALISE(GENERATED_MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GENERATED_MESHES_INITIALISE_GENERIC",ERR,ERROR)
    CALL EXITS("GENERATED_MESHES_INITIALISE_GENERIC")
    RETURN 1   
  END SUBROUTINE GENERATED_MESHES_INITIALISE_GENERIC
  
  !
  !================================================================================================================================
  !

  !>Intialises the generated meshes for an interface.
  SUBROUTINE GENERATED_MESHES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the generated meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%GENERATED_MESHES)) THEN
        CALL FLAG_ERROR("Interface generated meshes is already associated.",ERR,ERROR,*999)
      ELSE
        CALL GENERATED_MESHES_INITIALISE_GENERIC(INTERFACE%GENERATED_MESHES,ERR,ERROR,*999)
        INTERFACE%GENERATED_MESHES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESHES_INITIALISE_INTERFACE")
    RETURN
999 CALL ERRORS("GENERATED_MESHES_INITIALISE_INTERFACE",ERR,ERROR)
    CALL EXITS("GENERATED_MESHES_INITIALISE_INTERFACE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESHES_INITIALISE_INTERFACE
  
  !
  !================================================================================================================================
  !

  !>Intialises the generated meshes for a region.
  SUBROUTINE GENERATED_MESHES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%GENERATED_MESHES)) THEN
        CALL FLAG_ERROR("Region generated meshes is already associated.",ERR,ERROR,*999)
      ELSE
        CALL GENERATED_MESHES_INITIALISE_GENERIC(REGION%GENERATED_MESHES,ERR,ERROR,*999)
        REGION%GENERATED_MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
  
    CALL EXITS("GENERATED_MESHES_INITIALISE_REGION")
    RETURN
999 CALL ERRORS("GENERATED_MESHES_INITIALISE_REGION",ERR,ERROR)
    CALL EXITS("GENERATED_MESHES_INITIALISE_REGION")
    RETURN 1   
  END SUBROUTINE GENERATED_MESHES_INITIALISE_REGION
  
  !
  !================================================================================================================================
  ! 
  
  !>Updates the geometric field parameters from the initial nodal positions of the mesh. \see OPENCMISS::CMISSGeneratedMeshGeometricParametersCalculate
  SUBROUTINE GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE(FIELD,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the geometric parameters for
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<The mesh which is generated by the generated mesh \todo is this necessary???
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(ASSOCIATED(GENERATED_MESH)) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            CALL GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE(GENERATED_MESH%REGULAR_MESH,FIELD,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE(GENERATED_MESH%CYLINDER_MESH,FIELD,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented.",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE

  !
  !================================================================================================================================
  ! 

  !>Updates the geometric field parameters from the initial nodal positions of the regular mesh. Any derivative values for the nodes are calculated from an average straight line approximation.
  SUBROUTINE GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE(REGULAR_MESH,FIELD,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH !<A pointer to the regular mesh object
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! Local variables
    INTEGER(INTG) :: component_idx,global_np,global_np1,global_np2,nk,nk1,nk2,nl,nnl,np,np1,np2,ni,ny,&
      & DERIVATIVES_NUMBER_OF_LINES(8), TOTAL_NUMBER_OF_NODES_XI(3), node_idx(3),nd
    REAL(DP) :: DELTA(8),DELTA_COORD(3),MY_ORIGIN(3),MY_EXTENT(3),MESH_SIZE(3)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    !TODO: this subroutine can do with some more input parameter checking
    ! e.g. rectangular_cartesian coordinates checkings

    MY_ORIGIN=0.0_DP
    MY_EXTENT=0.0_DP
    MY_ORIGIN(1:REGULAR_MESH%MESH_DIMENSION)=REGULAR_MESH%ORIGIN
    MY_EXTENT(1:REGULAR_MESH%MESH_DIMENSION)=REGULAR_MESH%MAXIMUM_EXTENT
    MESH_SIZE=MY_EXTENT
    TOTAL_NUMBER_OF_NODES_XI=1
    DO ni=1,REGULAR_MESH%BASIS%NUMBER_OF_XI
      DELTA_COORD(ni)=MESH_SIZE(ni)/REAL(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni),DP)
      TOTAL_NUMBER_OF_NODES_XI(ni)=(REGULAR_MESH%BASIS%NUMBER_OF_NODES_XI(ni)-2)*REGULAR_MESH% &
        & NUMBER_OF_ELEMENTS_XI(ni)+REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)+1
    ENDDO !ni
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          FIELD_VARIABLE_COMPONENT=>FIELD_VARIABLE%COMPONENTS(component_idx)
          IF(FIELD_VARIABLE_COMPONENT%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
            ! update geometric parameters in this computational domain only
            DOMAIN=>FIELD_VARIABLE_COMPONENT%DOMAIN
            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
            DOMAIN_LINES=>DOMAIN%TOPOLOGY%LINES
            DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
              global_np=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
              node_idx(3)=(global_np-1)/(TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))+1
              node_idx(2)=MOD(global_np-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))/ &
                & TOTAL_NUMBER_OF_NODES_XI(1)+1
              node_idx(1)=MOD(MOD(global_np-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1)), &
                & TOTAL_NUMBER_OF_NODES_XI(1))+1
              ny=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,np)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ny, &
                & MY_ORIGIN(component_idx)+REAL(node_idx(component_idx)-1,DP)*DELTA_COORD(component_idx), &
                & ERR,ERROR,*999)
              ! calculate derivatives
              IF(DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES>1) THEN
                DERIVATIVES_NUMBER_OF_LINES=0
                DELTA=0.0_DP
                DO nnl=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES
                  nl=DOMAIN_NODES%NODES(np)%NODE_LINES(nnl)
                  np1=DOMAIN_LINES%LINES(nl)%NODES_IN_LINE(1)
                  global_np1=DOMAIN_NODES%NODES(np1)%GLOBAL_NUMBER
                  np2=DOMAIN_LINES%LINES(nl)%NODES_IN_LINE(DOMAIN_LINES%LINES(nl)%BASIS%NUMBER_OF_NODES)
                  global_np2=DOMAIN_NODES%NODES(np2)%GLOBAL_NUMBER
                  nk1=DOMAIN_LINES%LINES(nl)%DERIVATIVES_IN_LINE(2,1)
                  nk2=DOMAIN_LINES%LINES(nl)%DERIVATIVES_IN_LINE(2,DOMAIN_LINES%LINES(nl)%BASIS%NUMBER_OF_NODES)
                  !TODO: Adjust delta calculation for polar coordinate discontinuities
                  !TODO: this is hardcoded bit
                  SELECT CASE(component_idx)
                  CASE(1)
                    nd=MOD(global_np2-global_np1,TOTAL_NUMBER_OF_NODES_XI(component_idx))
                  CASE(2)
                    nd=MOD(global_np2-global_np1,TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)) &
                      & /TOTAL_NUMBER_OF_NODES_XI(1)
                  CASE(3)
                    nd=(global_np2-global_np1)/(TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2))
                  CASE DEFAULT
                    CALL FLAG_ERROR("Invalid component number.",ERR,ERROR,*999)
                  END SELECT
                  IF(np1==np) THEN
                    DERIVATIVES_NUMBER_OF_LINES(nk1)=DERIVATIVES_NUMBER_OF_LINES(nk1)+1
                    DELTA(nk1)=DELTA(nk1)+REAL(nd,DP)*DELTA_COORD(component_idx)
                  ELSE IF(np2==np) THEN
                    DERIVATIVES_NUMBER_OF_LINES(nk2)=DERIVATIVES_NUMBER_OF_LINES(nk2)+1
                    DELTA(nk2)=DELTA(nk2)
                    DELTA(nk2)=DELTA(nk2)+REAL(nd,DP)*DELTA_COORD(component_idx)
                  ELSE
                    !Error???
                  ENDIF
                ENDDO !nnl
                DO nk=1,8
                  IF(DERIVATIVES_NUMBER_OF_LINES(nk)>0) THEN
                    DELTA(nk)=DELTA(nk)/REAL(DERIVATIVES_NUMBER_OF_LINES(nk),DP)
                    ny=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(nk,np)
                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & ny,DELTA(nk),ERR,ERROR,*999)
                  ENDIF
                ENDDO !nk
              ENDIF
            ENDDO !np
          ELSE
            LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
              & " does not have node based interpolation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !component_idx
    !!TODO: do boundary nodes first then start the update to overlap computation and computation.
        CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
        CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
      ELSE
        LOCAL_ERROR="The standard field variable is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE

  !
  !================================================================================================================================
  ! 

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE(CYLINDER_MESH,FIELD,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH !<A pointer to the cylinder mesh object
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_TYPE),POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    INTEGER(INTG) :: NUMBER_ELEMENTS_XI(3),NUMBER_OF_NODES_XI(3)
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3),INTERPOLATION_TYPES(3)
    INTEGER(INTG) :: component_idx,xi_idx
    INTEGER(INTG) :: np,global_np,ny,nk
    INTEGER(INTG) :: NUMBER_OF_PLANAR_NODES,SCALING_TYPE
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: node_idx(3) ! holds r,theta,z indices
    REAL(DP) :: DELTA(3),DELTAi(3),POLAR_COORDS(3),RECT_COORDS(3)
    REAL(DP) :: CYLINDER_EXTENT(3),DERIV
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    ! calculate the total number of nodes in each xi direction
    IF(ASSOCIATED(CYLINDER_MESH%BASIS)) THEN
      BASIS=>CYLINDER_MESH%BASIS
      NUMBER_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
      NUMBER_OF_NODES_XI=BASIS%NUMBER_OF_NODES_XI
      DO xi_idx=1,3
        TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XI(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
      ENDDO
      TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
      ! calculate DELTAi now
      CYLINDER_EXTENT=CYLINDER_MESH%CYLINDER_EXTENT
      DELTA(1)=(CYLINDER_EXTENT(2)-CYLINDER_EXTENT(1))/NUMBER_ELEMENTS_XI(1)
      DELTA(2)=TWOPI/NUMBER_ELEMENTS_XI(2)
      DELTA(3)=CYLINDER_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
      DO xi_idx=1,3
        DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XI(xi_idx)-1)
      ENDDO
    ELSE
      CALL FLAG_ERROR("Cylinder mesh does not have a basis associated.",ERR,ERROR,*999)
    ENDIF

    ! assign to the field
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
          CALL FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,ERR,ERROR,*999)
          IF(SCALING_TYPE/=FIELD_UNIT_SCALING) &
            & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"  Note: If the cyinder looks wonky, set field scaling to&
            & unit scaling type.",ERR,ERROR,*999)
          NUMBER_OF_PLANAR_NODES=TOTAL_NUMBER_NODES_XI(1)*TOTAL_NUMBER_NODES_XI(2)
          DO component_idx=1,3
            INTERPOLATION_TYPES(component_idx)=FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE
          ENDDO
          IF(ALL(INTERPOLATION_TYPES==FIELD_NODE_BASED_INTERPOLATION)) THEN
            DOMAIN=>FIELD_VARIABLE%COMPONENTS(1)%DOMAIN ! just grab the first one
            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
            DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
              global_np=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
              ! calculate node_idx which will be used to calculate (r,theta,z) then (x,y,z)
              global_np=global_np-1 ! let's go 0-based index for a bit
              node_idx(3)=global_np/NUMBER_OF_PLANAR_NODES
              node_idx(2)=(global_np-(node_idx(3))*NUMBER_OF_PLANAR_NODES)/TOTAL_NUMBER_NODES_XI(1)
              node_idx(1)=MOD(global_np-(node_idx(3))*NUMBER_OF_PLANAR_NODES,TOTAL_NUMBER_NODES_XI(1))
              DO xi_idx=1,3
                POLAR_COORDS(xi_idx)=node_idx(xi_idx)*DELTAi(xi_idx)
              ENDDO
              POLAR_COORDS(1)=node_idx(1)*DELTAi(1)+CYLINDER_EXTENT(1) ! add the inner radius
              RECT_COORDS(1)=POLAR_COORDS(1)*COS(POLAR_COORDS(2))
              RECT_COORDS(2)=POLAR_COORDS(1)*SIN(POLAR_COORDS(2))
              RECT_COORDS(3)=POLAR_COORDS(3)
              RECT_COORDS=RECT_COORDS+CYLINDER_MESH%ORIGIN
              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                FIELD_VARIABLE_COMPONENT=>FIELD_VARIABLE%COMPONENTS(component_idx)
                ny=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,np)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ny, &
                  & RECT_COORDS(component_idx),ERR,ERROR,*999)
                ! Do derivatives: if there are derivatives, we can assume it's cubic hermite
                !   given that quadratic hermites are only used for collapsed hex elements,
                !   but NB mixed bases have to be handled (e.g. CH-CH-linear combinations)
                IF(DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES>1) THEN
                  ! Since I decided how xi 1,2,3 line up with the cylinder polar coordinates,
                  ! we know a priori that only some of the derivatives are nonzero (analytically).
                  ! NOTE: if hermite type used, should assign FIELD_UNIT_SCALING type for this to work
                  DO nk=2,FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES
                    SELECT CASE(DOMAIN_NODES%NODES(np)%GLOBAL_DERIVATIVE_INDEX(nk))
                    CASE(GLOBAL_DERIV_S1)
                      SELECT CASE(component_idx)
                      CASE(1)
                        DERIV=COS(POLAR_COORDS(2))*DELTA(1)
                      CASE(2)
                        DERIV=SIN(POLAR_COORDS(2))*DELTA(1)
                      CASE DEFAULT
                        DERIV=0.0_DP
                      END SELECT
                    CASE(GLOBAL_DERIV_S2)
                      SELECT CASE(component_idx)
                      CASE(1)
                        DERIV=-POLAR_COORDS(1)*SIN(POLAR_COORDS(2))*DELTA(2)
                      CASE(2)
                        DERIV=POLAR_COORDS(1)*COS(POLAR_COORDS(2))*DELTA(2)
                      CASE DEFAULT
                        DERIV=0.0_DP
                      END SELECT
                    CASE(GLOBAL_DERIV_S3)
                      IF(component_idx==3) THEN
                        DERIV=DELTA(3)
                      ELSE
                        DERIV=0.0_DP
                      ENDIF
                    CASE(GLOBAL_DERIV_S1_S2)
                      SELECT CASE(component_idx)
                      CASE(1)
                        DERIV=-SIN(POLAR_COORDS(2))*DELTA(1)*DELTA(2)
                      CASE(2)
                        DERIV=COS(POLAR_COORDS(2))*DELTA(1)*DELTA(2)
                      CASE DEFAULT
                        DERIV=0.0_DP
                      END SELECT
                    CASE DEFAULT  ! all other non-xy-planar cross derivatives
                      DERIV=0.0_DP
                    END SELECT
                    ! assign derivative
                    ny=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(nk,np)
                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                         & ny,DERIV,ERR,ERROR,*999)
                  ENDDO !nk
                ENDIF !derivatives
              ENDDO !component_idx
            ENDDO !np
          ELSE
            CALL FLAG_ERROR("All field variable components must have node-based interpolation.",ERR,ERROR,*999)
          ENDIF
          CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Geometric field must be three dimensional.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The standard field variable is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF

    ! all done
    IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)

    CALL EXITS("GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    CALL ERRORS("GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN 1  

  END SUBROUTINE GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE

  !
  !================================================================================================================================
  ! 

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_CYLINDER_SURFACES_GET(CYLINDER_MESH,SURFACE_TYPE,ELEMENTS,NODES,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH !<A pointer to the cylinder mesh object
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: ELEMENTS(:,:) !<On exit, contains the list of elements owning the surface (global numbers?)
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: NODES(:,:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: total_number_of_nodes,total_number_of_elements
    REAL(DP) :: delta(3),deltai(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GENERATED_MESH_CYLINDER_SURFACES_GET",ERR,ERROR,*999)

    ! let's go
    IF(ALLOCATED(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
      NUMBER_OF_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
      IF(ASSOCIATED(CYLINDER_MESH%BASIS)) THEN
        BASIS=>CYLINDER_MESH%BASIS
        IF(.NOT.ALLOCATED(ELEMENTS)) THEN
          IF(.NOT.ALLOCATED(NODES)) THEN
            NUMBER_OF_NODES_XI=BASIS%NUMBER_OF_NODES_XI
            ! build indices first (some of these are dummy arguments)
            CALL GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
              & cylinder_mesh%cylinder_extent,total_number_of_nodes,total_number_of_elements,NIDX,EIDX, &
              & delta,deltai,ERR,ERROR,*999)
            SELECT CASE(SURFACE_TYPE)
            CASE(GENERATED_MESH_CYLINDER_INNER_SURFACE)
              ALLOCATE(NODES(SIZE(NIDX,2),SIZE(NIDX,3)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              ALLOCATE(ELEMENTS(SIZE(EIDX,2),SIZE(EIDX,3)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENTS array.",ERR,ERROR,*999)
              NODES=NIDX(1,:,:)
              ELEMENTS=EIDX(1,:,:)
            CASE(GENERATED_MESH_CYLINDER_OUTER_SURFACE)
              ALLOCATE(NODES(SIZE(NIDX,2),SIZE(NIDX,3)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              ALLOCATE(ELEMENTS(SIZE(EIDX,2),SIZE(EIDX,3)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENTS array.",ERR,ERROR,*999)
              NODES=NIDX(SIZE(NIDX,1),:,:)
              ELEMENTS=EIDX(SIZE(EIDX,1),:,:)
            CASE(GENERATED_MESH_CYLINDER_TOP_SURFACE)
              ALLOCATE(NODES(SIZE(NIDX,1),SIZE(NIDX,2)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              ALLOCATE(ELEMENTS(SIZE(EIDX,1),SIZE(EIDX,2)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENTS array.",ERR,ERROR,*999)
              NODES=NIDX(:,:,SIZE(NIDX,3))
              ELEMENTS=EIDX(:,:,SIZE(EIDX,3))
            CASE(GENERATED_MESH_CYLINDER_BOTTOM_SURFACE)
              ALLOCATE(NODES(SIZE(NIDX,1),SIZE(NIDX,2)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              ALLOCATE(ELEMENTS(SIZE(EIDX,1),SIZE(EIDX,2)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENTS array.",ERR,ERROR,*999)
              NODES=NIDX(:,:,1)
              ELEMENTS=EIDX(:,:,1)
            CASE DEFAULT
              LOCAL_ERROR="The specified surface type of "//TRIM(NUMBER_TO_VSTRING(SURFACE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Output NODES array is already allocated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Output ELEMENTS array is already allocated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Cylinder mesh object does not have a basis associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Cylinder mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_CYLINDER_SURFACES_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CYLINDER_SURFACES_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_SURFACES_GET")
    RETURN 1  
  END SUBROUTINE GENERATED_MESH_CYLINDER_SURFACES_GET

  !
  !================================================================================================================================
  ! 

  !>Calculates the mesh topology information for a given cylinder (Not to be called by user)
  SUBROUTINE GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XI,CYLINDER_EXTENT, &
    & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAi,ERR,ERROR,*)
    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XI(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: CYLINDER_EXTENT(3)         !<inner & outer radii and height of cylinder
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_NODES    !<On exit, contains total number of nodes in cylinder mesh
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS !<On exit, contains total number of elements in cylinder mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: NIDX(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: EIDX(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: DELTA(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: DELTAi(3) !<Step sizes in each of (r,theta,z) for node (identical to DELTA if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) ! total number of nodes in each xi direction
    
    CALL ENTERS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES",ERR,ERROR,*999)
    
    ! Can skip most of the testing as this subroutine is only to be called by
    ! GENERATED_MESH_CYLINDER_CREATE_FINISH, which tests the input params.
    IF(.NOT.ALLOCATED(NIDX)) THEN
      IF(.NOT.ALLOCATED(EIDX)) THEN
        ! calculate DELTA and DELTAi
        DELTA(1)=(CYLINDER_EXTENT(2)-CYLINDER_EXTENT(1))/NUMBER_ELEMENTS_XI(1)
        DELTA(2)=TWOPI/NUMBER_ELEMENTS_XI(2)
        DELTA(3)=CYLINDER_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
        DO xi_idx=1,3
          DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XI(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XI(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
        TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)
        
        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NIDX array.",ERR,ERROR,*999)
        NN=0
        DO nn3=1,TOTAL_NUMBER_NODES_XI(3)
          DO nn2=1,TOTAL_NUMBER_NODES_XI(2)
            DO nn1=1,TOTAL_NUMBER_NODES_XI(1)
              NN=NN+1
              NIDX(nn1,nn2,nn3)=NN
            ENDDO ! nn1
          ENDDO ! nn2
        ENDDO ! nn3
        TOTAL_NUMBER_OF_NODES=NN

        ! now do EIDX
        ALLOCATE(EIDX(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2),NUMBER_ELEMENTS_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate EIDX array.",ERR,ERROR,*999)
        NE=0
        DO ne3=1,NUMBER_ELEMENTS_XI(3)
          DO ne2=1,NUMBER_ELEMENTS_XI(2)
            DO ne1=1,NUMBER_ELEMENTS_XI(1)
              NE=NE+1
              EIDX(ne1,ne2,ne3)=NE
            ENDDO
          ENDDO
        ENDDO
        TOTAL_NUMBER_OF_ELEMENTS=NE
      
      ELSE
        CALL FLAG_ERROR("NIDX array is already allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("EIDX array is already allocated.",ERR,error,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES

END MODULE GENERATED_MESH_ROUTINES
