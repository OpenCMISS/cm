!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all generated mesh routines.
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
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS
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
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_MESH_TYPE=5 !<An ellipsoid generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  !>@}

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces GENERATED_MESH_ROUTINES::GeneratedMeshCylinderSurfaces
  !> \brief Generated mesh cylinder type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_INNER_SURFACE=1  !<Inner surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_OUTER_SURFACE=2  !<Outer surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_TOP_SURFACE=3    !<Top surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_BOTTOM_SURFACE=4 !<Bottom surface of the cylinder. \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
  !>@}

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshEllipsoidSurfaces GENERATED_MESH_ROUTINES::GeneratedMeshEllipsoidSurfaces
  !> \brief Generated mesh ellipsoid type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_INNER_SURFACE=5  !<Inner surface of the ellipsoid. \see GENERATED_MESH_ROUTINES_GeneratedMeshEllipsoidSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_OUTER_SURFACE=6  !<Outer surface of the ellipsoid. \see GENERATED_MESH_ROUTINES_GeneratedMeshEllipsoidSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_TOP_SURFACE=7    !<Top surface of the ellipsoid. \see GENERATED_MESH_ROUTINES_GeneratedMeshEllipsoidSurfaces,GENERATED_MESH_ROUTINES
  !>@}

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces GENERATED_MESH_ROUTINES::GeneratedMeshRegularSurfaces
  !> \brief Generated mesh regular type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_LEFT_SURFACE=8    !<Left surface of the regular mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_RIGHT_SURFACE=9   !<Right surface of the regular mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_TOP_SURFACE=10    !<Top surface of the regular mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BOTTOM_SURFACE=1  !<Bottom surface of the regular mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_FRONT_SURFACE=12  !<Front surface of the regular mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BACK_SURFACE=13   !<Back surface of the regular mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
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
  PUBLIC GENERATED_MESH_ELLIPSOID_MESH_TYPE
  PUBLIC GENERATED_MESH_CYLINDER_INNER_SURFACE,GENERATED_MESH_CYLINDER_OUTER_SURFACE
  PUBLIC GENERATED_MESH_CYLINDER_TOP_SURFACE,GENERATED_MESH_CYLINDER_BOTTOM_SURFACE
  PUBLIC GENERATED_MESH_ELLIPSOID_INNER_SURFACE,GENERATED_MESH_ELLIPSOID_OUTER_SURFACE
  PUBLIC GENERATED_MESH_ELLIPSOID_TOP_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_LEFT_SURFACE,GENERATED_MESH_REGULAR_RIGHT_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_TOP_SURFACE,GENERATED_MESH_REGULAR_BOTTOM_SURFACE
  PUBLIC GENERATED_MESH_REGULAR_FRONT_SURFACE,GENERATED_MESH_REGULAR_BACK_SURFACE
  PUBLIC GENERATED_MESHES_INITIALISE,GENERATED_MESHES_FINALISE

  PUBLIC GENERATED_MESH_BASE_VECTORS_SET

  PUBLIC GENERATED_MESH_COORDINATE_SYSTEM_GET

  PUBLIC GENERATED_MESH_CREATE_START,GENERATED_MESH_CREATE_FINISH

  PUBLIC GENERATED_MESH_DESTROY

  PUBLIC GENERATED_MESH_BASIS_SET,GENERATED_MESH_EXTENT_SET,GENERATED_MESH_NUMBER_OF_ELEMENTS_SET,GENERATED_MESH_ORIGIN_SET, &
    & GENERATED_MESH_TYPE_SET, GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE

  PUBLIC GENERATED_MESH_BASIS_GET,GENERATED_MESH_EXTENT_GET,GENERATED_MESH_NUMBER_OF_ELEMENTS_GET,GENERATED_MESH_ORIGIN_GET,&
    & GENERATED_MESH_TYPE_GET

  PUBLIC GENERATED_MESH_REGION_GET

  PUBLIC GENERATED_MESH_USER_NUMBER_FIND
  PUBLIC GENERATED_MESH_SURFACE_GET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the basis of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBasisGet
  SUBROUTINE GENERATED_MESH_BASIS_GET(GENERATED_MESH,BASES,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the bases of
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:) !<On return, the bases of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: basis_idx,NUM_BASES

    CALL ENTERS("GENERATED_MESH_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) THEN
              NUM_BASES=SIZE(GENERATED_MESH%REGULAR_MESH%BASES)
              ALLOCATE(BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASES(basis_idx)%PTR=>GENERATED_MESH%REGULAR_MESH%BASES(basis_idx)%PTR
              ENDDO
            ELSE
              CALL FLAG_ERROR("Generated mesh bases are not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Generated mesh regular mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) THEN
              NUM_BASES=SIZE(GENERATED_MESH%CYLINDER_MESH%BASES)
              ALLOCATE(BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASES(basis_idx)%PTR=>GENERATED_MESH%CYLINDER_MESH%BASES(basis_idx)%PTR
              ENDDO
            ELSE
              CALL FLAG_ERROR("Generated mesh bases are not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Generated mesh cylinder mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) THEN
              NUM_BASES=SIZE(GENERATED_MESH%ELLIPSOID_MESH%BASES)
              ALLOCATE(BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASES(basis_idx)%PTR=>GENERATED_MESH%ELLIPSOID_MESH%BASES(basis_idx)%PTR
              ENDDO
            ELSE
              CALL FLAG_ERROR("Generated mesh bases are not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Generated mesh ellipsoid mesh is not associated.",ERR,ERROR,*999)
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
  SUBROUTINE GENERATED_MESH_BASIS_SET(GENERATED_MESH,BASES,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the basis of
    TYPE(BASIS_PTR_TYPE) :: BASES(:) !<An array of pointers to the basis to generate the mesh with
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: COORDINATE_DIMENSION,basis_idx, NUM_BASES, NUM_XI, BASIS_TYPE
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,COORDINATE_DIMENSION,ERR,ERROR,*999)
        NUM_BASES=SIZE(BASES)
        NUM_XI=BASES(1)%PTR%NUMBER_OF_XI
        BASIS_TYPE=BASES(1)%PTR%TYPE
        DO basis_idx=2,NUM_BASES
          IF(BASES(basis_idx)%PTR%NUMBER_OF_XI /= NUM_XI) THEN
            CALL FLAG_ERROR("All bases must have the same number of xi.",ERR,ERROR,*999)
          ENDIF
          IF(BASES(basis_idx)%PTR%TYPE /= BASIS_TYPE) THEN
            CALL FLAG_ERROR("Using different basis types is not supported for generated meshes.",ERR,ERROR,*999)
          ENDIF
        ENDDO
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS)) THEN
              CALL FLAG_ERROR("Can not reset the basis if base vectors have been specified.",ERR,ERROR,*999)
            ELSE
              IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASES)) DEALLOCATE(GENERATED_MESH%REGULAR_MESH%BASES)
              ALLOCATE(GENERATED_MESH%REGULAR_MESH%BASES(NUM_BASES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                IF(ASSOCIATED(BASES(basis_idx)%PTR)) THEN
                  IF(BASES(basis_idx)%PTR%NUMBER_OF_XI<=COORDINATE_DIMENSION) THEN
                    GENERATED_MESH%REGULAR_MESH%BASES(basis_idx)%PTR=>BASES(basis_idx)%PTR
                  ELSE
                    LOCAL_ERROR="The basis number of xi dimensions of "// &
                      & TRIM(NUMBER_TO_VSTRING(BASES(basis_idx)%PTR%NUMBER_OF_XI,"*",ERR,ERROR))// &
                      & " is invalid. The number of xi dimensions must be <= the number of coordinate dimensions of "// &
                      & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The basis with index "//TRIM(NUMBER_TO_VSTRING(basis_idx,"*",ERR,ERROR))// &
                    & " is not associated."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO
            ENDIF
          ELSE
            CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%CYLINDER_MESH%BASES)) DEALLOCATE(GENERATED_MESH%CYLINDER_MESH%BASES)
            ALLOCATE(GENERATED_MESH%CYLINDER_MESH%BASES(NUM_BASES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
            DO basis_idx=1,NUM_BASES
              IF(ASSOCIATED(BASES(basis_idx)%PTR)) THEN
                GENERATED_MESH%CYLINDER_MESH%BASES(basis_idx)%PTR=>BASES(basis_idx)%PTR
              ELSE
                LOCAL_ERROR="The basis with index "//TRIM(NUMBER_TO_VSTRING(basis_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            CALL FLAG_ERROR("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
            IF(ALLOCATED(GENERATED_MESH%ELLIPSOID_MESH%BASES)) DEALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%BASES)
            ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%BASES(NUM_BASES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
            DO basis_idx=1,NUM_BASES
              IF(ASSOCIATED(BASES(basis_idx)%PTR)) THEN
                GENERATED_MESH%ELLIPSOID_MESH%BASES(basis_idx)%PTR=>BASES(basis_idx)%PTR
              ELSE
                LOCAL_ERROR="The basis with index "//TRIM(NUMBER_TO_VSTRING(basis_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            CALL FLAG_ERROR("Ellpsoid generated mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
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

  !>Sets/changes the base vectors of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshBaseVectorsSet
  SUBROUTINE GENERATED_MESH_BASE_VECTORS_SET(GENERATED_MESH,BASE_VECTORS,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the base vectors fo
    REAL(DP), INTENT(IN) :: BASE_VECTORS(:,:) !<BASE_VECTORS(coordinate_idx,xi_idx). The base vectors for the generated mesh to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: COORDINATE_DIMENSION
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_BASE_VECTORS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,COORDINATE_DIMENSION,ERR,ERROR,*999)
        IF(SIZE(BASE_VECTORS,1)==COORDINATE_DIMENSION) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
              BASES=>GENERATED_MESH%REGULAR_MESH%BASES
              IF(ASSOCIATED(BASES)) THEN
                BASIS=>BASES(1)%PTR !Bases should all have same number of xi
                IF(ASSOCIATED(BASIS)) THEN
                  IF(SIZE(BASE_VECTORS,2)==BASIS%NUMBER_OF_XI) THEN
                    IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS)) DEALLOCATE(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS)
                    ALLOCATE(GENERATED_MESH%REGULAR_MESH%BASE_VECTORS(SIZE(BASE_VECTORS,1),SIZE(BASE_VECTORS,2)),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate base vectors.",ERR,ERROR,*999)
                    GENERATED_MESH%REGULAR_MESH%BASE_VECTORS=BASE_VECTORS
                  ELSE
                    LOCAL_ERROR="The size of the second dimension of base vectors of "// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(BASE_VECTORS,2),"*",ERR,ERROR))// &
                      & " is invalid. The second dimension size must match the number of mesh dimensions of "// &
                      & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Bases are not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("You must set the generated mesh basis before setting base vectors.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          LOCAL_ERROR="The size of the first dimension of base vectors of "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(BASE_VECTORS,1),"*",ERR,ERROR))// &
            & " is invalid. The first dimension size must match the coordinate system dimension of "// &
            & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_BASE_VECTORS_SET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_BASE_VECTORS_SET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_BASE_VECTORS_SET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_BASE_VECTORS_SET

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for a generated mesh accounting for regions and interfaces
  SUBROUTINE GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On return, the generated meshes coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(REGION_TYPE), POINTER :: REGION,PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
        CALL FLAG_ERROR("Coordinate system is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          COORDINATE_SYSTEM=>REGION%COORDINATE_SYSTEM
          IF(.NOT.ASSOCIATED(COORDINATE_SYSTEM)) THEN
            LOCAL_ERROR="The coordinate system is not associated for generated mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of region number "// &
              & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          NULLIFY(INTERFACE)
          INTERFACE=>GENERATED_MESH%INTERFACE
          IF(ASSOCIATED(INTERFACE)) THEN
            PARENT_REGION=>INTERFACE%PARENT_REGION
            IF(ASSOCIATED(PARENT_REGION)) THEN
              COORDINATE_SYSTEM=>PARENT_REGION%COORDINATE_SYSTEM
              IF(.NOT.ASSOCIATED(COORDINATE_SYSTEM)) THEN
                LOCAL_ERROR="The coordinate system is not associated for generated mesh number "// &
                  & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of interface number "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" of parent region number "// &
                  & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The parent region not associated for generated mesh number "// &
                & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of interface number "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The region or interface is not associated for generated mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_COORDINATE_SYSTEM_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_COORDINATE_SYSTEM_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_COORDINATE_SYSTEM_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_COORDINATE_SYSTEM_GET

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
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
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
          GENERATED_MESH%GENERATED_MESH_FINISHED=.TRUE.
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
            GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES-1
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
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        EXTENT=GENERATED_MESH%ELLIPSOID_MESH%ELLIPSOID_EXTENT
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
    INTEGER(INTG) :: COORDINATE_DIMENSION
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_EXTENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,COORDINATE_DIMENSION,ERR,ERROR,*999)
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(SIZE(EXTENT,1)==COORDINATE_DIMENSION) THEN
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
              IF(L2NORM(EXTENT)>ZERO_TOLERANCE) THEN
                IF(ALLOCATED(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)) &
                    & DEALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)
                ALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT(COORDINATE_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate maximum extent.",ERR,ERROR,*999)
                GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT=EXTENT
              ELSE
                CALL FLAG_ERROR("The norm of the mesh extent is zero.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Regular generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The extent size of "//TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))// &
                & " is invalid. The extent size must match the coordinate system dimension of "// &
                & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_POLAR_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
          IF(SIZE(EXTENT,1)==COORDINATE_DIMENSION) THEN
            IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
              ALLOCATE(GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT(SIZE(EXTENT)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate maximum extent.",ERR,ERROR,*999)
              GENERATED_MESH%CYLINDER_MESH%CYLINDER_EXTENT=EXTENT
            ELSE
              CALL FLAG_ERROR("Cylinder generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The extent size of "//TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))// &
                & " is invalid. The extent size must match the coordinate system dimension of "// &
                & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF((SIZE(EXTENT,1)-1)==COORDINATE_DIMENSION) THEN
            IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
              ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%ELLIPSOID_EXTENT(SIZE(EXTENT)),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate maximum extent.",ERR,ERROR,*999)
              GENERATED_MESH%ELLIPSOID_MESH%ELLIPSOID_EXTENT=EXTENT
            ELSE
              CALL FLAG_ERROR("Ellipsoid generated mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The extent size of "//TRIM(NUMBER_TO_VSTRING(SIZE(EXTENT,1),"*",ERR,ERROR))// &
                & " is invalid. The extent size must be equal one plus the coordinate system dimension of "// &
                & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The generated mesh mesh type of "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
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

  !>Get one of the surfaces of a generated mesh. \see OPENCMISS::CMISSGeneratedMeshSurfaceGet
  SUBROUTINE GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the type of
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<The surface you are interested in
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES (:) !<The nodes on the specified surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<The normal outward pointing xi direction
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR
!     INTEGER(INTG), ALLOCATABLE :: NODES(:)


    CALL ENTERS("GENERATED_MESH_SURFACE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH
        CALL GENERATED_MESH_REGULAR_SURFACE_GET(REGULAR_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,ERR,ERROR,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        CYLINDER_MESH=>GENERATED_MESH%CYLINDER_MESH
        CALL GENERATED_MESH_CYLINDER_SURFACE_GET(CYLINDER_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,ERR,ERROR,*999)
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        ELLIPSOID_MESH=>GENERATED_MESH%ELLIPSOID_MESH
        CALL GENERATED_MESH_ELLIPSOID_SURFACE_GET(ELLIPSOID_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI, &
            & ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_SURFACE_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_SURFACE_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_SURFACE_GET")
    RETURN
  END SUBROUTINE GENERATED_MESH_SURFACE_GET

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
      CALL GENERATED_MESH_ELLIPSOID_FINALISE(GENERATED_MESH%ELLIPSOID_MESH,ERR,ERROR,*999)
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
      NULLIFY(GENERATED_MESH%ELLIPSOID_MESH)
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
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(SIZE(NUMBER_OF_ELEMENTS,1)>=SIZE(GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI,1)) THEN
          NUMBER_OF_ELEMENTS=GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
          LOCAL_ERROR="The size of NUMBER_OF_ELEMENTS is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))//"."
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
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH
          IF(ASSOCIATED(REGULAR_MESH)) THEN
            BASIS=>REGULAR_MESH%BASES(1)%PTR !Number of xi will be the same for all bases
            IF(ASSOCIATED(BASIS)) THEN
              IF(SIZE(NUMBER_OF_ELEMENTS_XI,1)==BASIS%NUMBER_OF_XI) THEN
                IF(ALL(NUMBER_OF_ELEMENTS_XI>0)) THEN
                  IF(ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) DEALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)
                  ALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(BASIS%NUMBER_OF_XI),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                  REGULAR_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI
                ELSE
                  CALL FLAG_ERROR("Must have 1 or more elements in all directions.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The number of elements xi size of "// &
                  & TRIM(NUMBER_TO_VSTRING(SIZE(NUMBER_OF_ELEMENTS_XI,1),"*",ERR,ERROR))// &
                  & " is invalid. The number of elements xi size must match the basis number of xi dimensions of "// &
                  & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Must set the generated mesh basis before setting the number of elements.",ERR,ERROR,*999)
            ENDIF
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
        CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
            ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI(SIZE(NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
            GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI
          ELSE
            CALL FLAG_ERROR("Ellipsoid generated mesh is not associated.",ERR,ERROR,*999)
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
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(SIZE(ORIGIN,1)>=SIZE(GENERATED_MESH%ELLIPSOID_MESH%ORIGIN,1)) THEN
          ORIGIN=GENERATED_MESH%ELLIPSOID_MESH%ORIGIN
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
    INTEGER(INTG) :: COORDINATE_DIMENSION
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_ORIGIN_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished.",ERR,ERROR,*999)
      ELSE
        NULLIFY(COORDINATE_SYSTEM)
        CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
        CALL COORDINATE_SYSTEM_DIMENSION_GET(COORDINATE_SYSTEM,COORDINATE_DIMENSION,ERR,ERROR,*999)
        IF(SIZE(ORIGIN,1)==COORDINATE_DIMENSION) THEN
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
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
              IF(SIZE(ORIGIN,1)==3) THEN
                IF(.NOT.ALLOCATED(GENERATED_MESH%ELLIPSOID_MESH%ORIGIN)) THEN
                  ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH%ORIGIN(SIZE(ORIGIN)),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Ellipsoid generated mesh is only supported for 3D.",ERR,ERROR,*999)
              ENDIF
              GENERATED_MESH%ELLIPSOID_MESH%ORIGIN=ORIGIN
            ELSE
              CALL FLAG_ERROR("Ellipsoid generated mesh is not associated.",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="The generated mesh mesh type of "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          LOCAL_ERROR="The origin size of "//TRIM(NUMBER_TO_VSTRING(SIZE(ORIGIN,1),"*",ERR,ERROR))// &
            & " is invalid. The extent size must match the coordinate system dimension of "// &
            & TRIM(NUMBER_TO_VSTRING(COORDINATE_DIMENSION,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
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
    INTEGER(INTG) :: coordinate_idx,COUNT,ELEMENT_FACTOR,grid_ne,GRID_NUMBER_OF_ELEMENTS,ni,ne,ne1,ne2,ne3,nn,nn1,nn2,nn3,np, &
      & NUMBER_OF_ELEMENTS_XI(3),TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_OF_NODES,NUMBER_OF_CORNER_NODES, &
      & TOTAL_NUMBER_OF_ELEMENTS,xi_idx,NUM_BASES,basis_idx,BASIS_NUMBER_OF_NODES
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:),ELEMENT_NODES_USER_NUMBERS(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      NULLIFY(COORDINATE_SYSTEM)
      CALL GENERATED_MESH_COORDINATE_SYSTEM_GET(GENERATED_MESH,COORDINATE_SYSTEM,ERR,ERROR,*999)
      REGION=>GENERATED_MESH%REGION
      INTERFACE=>GENERATED_MESH%INTERFACE
      REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH
      IF(ASSOCIATED(REGULAR_MESH)) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(ALLOCATED(REGULAR_MESH%BASES)) THEN
            !Use first basis to get number of xi
            BASIS=>REGULAR_MESH%BASES(1)%PTR
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE)
              IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED)) &
                & CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
              !Determine the coordinate system and create the regular mesh for that system
              REGULAR_MESH%COORDINATE_DIMENSION=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              REGULAR_MESH%MESH_DIMENSION=BASIS%NUMBER_OF_XI
              IF(.NOT.ALLOCATED(REGULAR_MESH%ORIGIN)) THEN
                ALLOCATE(REGULAR_MESH%ORIGIN(REGULAR_MESH%COORDINATE_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
                REGULAR_MESH%ORIGIN=0.0_DP
              ENDIF
              IF(.NOT.ALLOCATED(REGULAR_MESH%MAXIMUM_EXTENT)) THEN
                ALLOCATE(REGULAR_MESH%MAXIMUM_EXTENT(REGULAR_MESH%COORDINATE_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate maximum extent.",ERR,ERROR,*999)
                REGULAR_MESH%MAXIMUM_EXTENT=1.0_DP
              ENDIF
              IF(.NOT.ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
                ALLOCATE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(REGULAR_MESH%MESH_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                REGULAR_MESH%NUMBER_OF_ELEMENTS_XI=1
              ENDIF
              IF(ALLOCATED(REGULAR_MESH%BASE_VECTORS)) THEN
                !!TODO: Check base vectors
              ELSE
                !Calculate base vectors
                ALLOCATE(REGULAR_MESH%BASE_VECTORS(REGULAR_MESH%COORDINATE_DIMENSION,REGULAR_MESH%MESH_DIMENSION),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                REGULAR_MESH%BASE_VECTORS=0.0_DP
                IF(REGULAR_MESH%MESH_DIMENSION==1) THEN
                  !The base vector is just the extent vector
                  REGULAR_MESH%BASE_VECTORS(:,1)=REGULAR_MESH%MAXIMUM_EXTENT
                ELSE
                  IF(REGULAR_MESH%MESH_DIMENSION<REGULAR_MESH%COORDINATE_DIMENSION) THEN
                    !Find the first number of mesh dimensions for which the extent is non-zero.
                    COUNT=0
                    coordinate_idx=1
                    DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                      DO WHILE(coordinate_idx<=REGULAR_MESH%COORDINATE_DIMENSION)
                        IF(ABS(REGULAR_MESH%MAXIMUM_EXTENT(coordinate_idx))>ZERO_TOLERANCE) THEN
                          REGULAR_MESH%BASE_VECTORS(coordinate_idx,xi_idx)=REGULAR_MESH%MAXIMUM_EXTENT(xi_idx)
                          COUNT=COUNT+1
                        ELSE
                          coordinate_idx=coordinate_idx+1
                        ENDIF
                      ENDDO
                    ENDDO !xi_idx
                    IF(COUNT/=REGULAR_MESH%MESH_DIMENSION)  &
                      & CALL FLAG_ERROR("Invalid mesh extent. There number of non-zero components is < the mesh dimension.", &
                      & ERR,ERROR,*999)
                  ELSE IF(REGULAR_MESH%MESH_DIMENSION==REGULAR_MESH%COORDINATE_DIMENSION) THEN
                    !The default base vectors are aligned with the coordinate vectors
                    DO coordinate_idx=1,REGULAR_MESH%COORDINATE_DIMENSION
                      REGULAR_MESH%BASE_VECTORS(coordinate_idx,coordinate_idx)=REGULAR_MESH%MAXIMUM_EXTENT(coordinate_idx)
                    ENDDO !coordinate_idx
                  ELSE
                    CALL FLAG_ERROR("The mesh dimension is greater than the coordinate dimension.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDIF
              !Calculate the sizes of a regular grid of elements with the appropriate number of basis nodes in each dimension of
              !the grid element
              TOTAL_NUMBER_OF_NODES=1
              GRID_NUMBER_OF_ELEMENTS=1
              NUMBER_OF_ELEMENTS_XI=1
              NUM_BASES=SIZE(REGULAR_MESH%BASES)
              DO ni=1,REGULAR_MESH%MESH_DIMENSION
                !Set total number of nodes to corner nodes only
                TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES*(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)+1)
                NUMBER_OF_ELEMENTS_XI(ni)=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)
                GRID_NUMBER_OF_ELEMENTS=GRID_NUMBER_OF_ELEMENTS*REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)
              ENDDO
              NUMBER_OF_CORNER_NODES=TOTAL_NUMBER_OF_NODES
              !Add extra nodes for each basis
              !Will end up with some duplicate nodes if bases have the same interpolation in one direction
              DO basis_idx=1,NUM_BASES
                BASIS=>REGULAR_MESH%BASES(basis_idx)%PTR
                BASIS_NUMBER_OF_NODES=1
                DO ni=1,REGULAR_MESH%MESH_DIMENSION
                  BASIS_NUMBER_OF_NODES=BASIS_NUMBER_OF_NODES*((BASIS%NUMBER_OF_NODES_XIC(ni)-1)* &
                      & REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)+1)
                ENDDO
                TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES-NUMBER_OF_CORNER_NODES
              ENDDO
              !Compute the element factor i.e., the number of sub elements each grid element will be split into.
              IF(BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                ELEMENT_FACTOR=1
              ELSE
                SELECT CASE(REGULAR_MESH%MESH_DIMENSION)
                CASE(1)
                  ELEMENT_FACTOR=1
                CASE(2)
                  ELEMENT_FACTOR=2
                CASE(3)
                  ELEMENT_FACTOR=6
                CASE DEFAULT
                  LOCAL_ERROR="The mesh dimension dimension of "// &
                    & TRIM(NUMBER_TO_VSTRING(REGULAR_MESH%MESH_DIMENSION,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
              TOTAL_NUMBER_OF_ELEMENTS=ELEMENT_FACTOR*GRID_NUMBER_OF_ELEMENTS
              !Create the default node set
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
                CALL MESH_CREATE_START(MESH_USER_NUMBER,REGION,REGULAR_MESH%MESH_DIMENSION,GENERATED_MESH%MESH, &
                  & ERR,ERROR,*999)
              ELSE
                CALL MESH_CREATE_START(MESH_USER_NUMBER,INTERFACE,REGULAR_MESH%MESH_DIMENSION,GENERATED_MESH%MESH, &
                  & ERR,ERROR,*999)
              ENDIF
              !Set the number of mesh components
              CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,NUM_BASES,ERR,ERROR,*999)
              !Create the elements
              CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
              DO basis_idx=1,NUM_BASES
                BASIS=>REGULAR_MESH%BASES(basis_idx)%PTR
                !Get number of nodes in each xi direction for this basis
                DO ni=1,REGULAR_MESH%MESH_DIMENSION
                  TOTAL_NUMBER_OF_NODES_XI(ni)=(BASIS%NUMBER_OF_NODES_XIC(ni)-1)*REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(ni)+1
                ENDDO
                NULLIFY(MESH_ELEMENTS)
                CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,basis_idx,BASIS,MESH_ELEMENTS,ERR,ERROR,*999)
                !Set the elements for the regular mesh
                IF (ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
                ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF (ALLOCATED(ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(ELEMENT_NODES_USER_NUMBERS)
                ALLOCATE(ELEMENT_NODES_USER_NUMBERS(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element nodes.",ERR,ERROR,*999)
                !Step in the xi(3) direction
                DO ne3=1,NUMBER_OF_ELEMENTS_XI(3)+1
                  DO ne2=1,NUMBER_OF_ELEMENTS_XI(2)+1
                    DO ne1=1,NUMBER_OF_ELEMENTS_XI(1)+1
                      IF(BASIS%NUMBER_OF_XI<3.OR.ne3<=NUMBER_OF_ELEMENTS_XI(3)) THEN
                        IF(BASIS%NUMBER_OF_XI<2.OR.ne2<=NUMBER_OF_ELEMENTS_XI(2)) THEN
                          IF(ne1<=NUMBER_OF_ELEMENTS_XI(1)) THEN
                            grid_ne=ne1
                            np=1+(ne1-1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)
                            IF(BASIS%NUMBER_OF_XI>1) THEN
                              grid_ne=grid_ne+(ne2-1)*NUMBER_OF_ELEMENTS_XI(1)
                              np=np+(ne2-1)*TOTAL_NUMBER_OF_NODES_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)
                              IF(BASIS%NUMBER_OF_XI>2) THEN
                                grid_ne=grid_ne+(ne3-1)*NUMBER_OF_ELEMENTS_XI(1)*NUMBER_OF_ELEMENTS_XI(2)
                                np=np+(ne3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)* &
                                  & (BASIS%NUMBER_OF_NODES_XIC(3)-1)
                              ENDIF
                            ENDIF
                            IF(BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
                              !Lagrange Hermite TP elements
                              ne=grid_ne
                              nn=0
                              DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
                                nn=nn+1
                                ELEMENT_NODES(nn)=np+(nn1-1)
                              ENDDO !nn1
                              IF(BASIS%NUMBER_OF_XI>1) THEN
                                DO nn2=2,BASIS%NUMBER_OF_NODES_XIC(2)
                                  DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
                                    nn=nn+1
                                    ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ENDDO !nn1
                                ENDDO !nn2
                                IF(BASIS%NUMBER_OF_XI>2) THEN
                                  DO nn3=2,BASIS%NUMBER_OF_NODES_XIC(3)
                                    DO nn2=1,BASIS%NUMBER_OF_NODES_XIC(2)
                                      DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
                                        nn=nn+1
                                        ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)+ &
                                            & (nn3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                      ENDDO !nn1
                                    ENDDO !nn2
                                  ENDDO !nn3
                                ENDIF
                              ENDIF
                              CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                  & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                              CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                  & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                            ELSE
                              !Simplex elements
                              SELECT CASE(BASIS%NUMBER_OF_XI)
                              CASE(1)
                                !Line element
                                ne=grid_ne
                                nn=0
                                DO nn1=1,BASIS%NUMBER_OF_NODES_XIC(1)
                                  nn=nn+1
                                  ELEMENT_NODES(nn)=np+(nn1-1)
                                ENDDO !nn1
                                CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                    & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                    & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                              CASE(2)
                                !Triangular element
                                !Break the grid square element into 2 triangles. The 2 triangles are
                                !Element 1: vertices {(0,0);(1,0);(1,1)}
                                !Element 2: vertices {(0,0);(0,1);(1,1)}
                                SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
                                CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
                                  !First sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+1
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+1
                                  ELEMENT_NODES(3)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Second sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+2
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+1
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2
                                  ELEMENT_NODES(3)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+1
                                  ELEMENT_NODES(5)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Second sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+2
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(5)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                               CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+1
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3
                                  ELEMENT_NODES(3)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+1
                                  ELEMENT_NODES(5)=np+2
                                  ELEMENT_NODES(6)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(7)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(8)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(10)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Second sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+2
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(5)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(7)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(8)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(10)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                CASE DEFAULT
                                  LOCAL_ERROR="The simplex basis interpolation order of "// &
                                    & TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(1),"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                              CASE(3)
                                !Tetrahedra element
                                !Break the grid cube element into 6 tetrahedra (so that we have a break down the main diagonal of the
                                !cube in order to allow for the middle node in quadratics to be included). The 6 tetrahedra are
                                !Element 1: vertices {(0,0,0);(1,0,0);(1,1,0);(1,1,1)}
                                !Element 2: vertices {(0,0,0);(0,1,0);(1,1,0);(1,1,1)}
                                !Element 3: vertices {(0,0,0);(1,0,0);(1,0,1);(1,1,1)}
                                !Element 4: vertices {(0,0,0);(0,0,1);(1,0,1);(1,1,1)}
                                !Element 5: vertices {(0,0,0);(0,1,0);(0,1,1);(1,1,1)}
                                !Element 6: vertices {(0,0,0);(0,0,1);(0,1,1);(1,1,1)}
                                SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
                                CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
                                  !First sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+1
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+1
                                  ELEMENT_NODES(3)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Second sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+2
                                   ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Third sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+3
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+1
                                  ELEMENT_NODES(3)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Fourth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+4
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(3)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Fifth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+5
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Sixth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+6
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(3)=np+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+1
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2
                                  ELEMENT_NODES(3)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+1
                                  ELEMENT_NODES(6)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(9)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Second sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+2
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(9)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Third sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+3
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2
                                  ELEMENT_NODES(3)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+1
                                  ELEMENT_NODES(6)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Fourth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+4
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(3)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(6)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Fifth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+5
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Sixth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+6
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(3)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(6)=np+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
                                  !First sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+1
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3
                                  ELEMENT_NODES(3)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+1
                                  ELEMENT_NODES(6)=np+2
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(8)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(11)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(12)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(13)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(14)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(15)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(16)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(17)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(18)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(19)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(20)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Second sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+2
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(4)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(8)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(11)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(12)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(13)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(14)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(15)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(16)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(17)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(18)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(19)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(20)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Third sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+3
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3
                                  ELEMENT_NODES(3)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+1
                                  ELEMENT_NODES(6)=np+2
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(11)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(12)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(13)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(14)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(15)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(16)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(17)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(18)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(19)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(20)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Fourth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+4
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(3)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(6)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(7)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(11)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(12)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(13)=np+3+TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(14)=np+3+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(15)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(16)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(17)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(18)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(19)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(20)=np+2+TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Fifth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+5
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(3)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(6)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)
                                  ELEMENT_NODES(7)=np+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(11)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(12)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(13)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(14)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(15)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(16)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(17)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(18)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(19)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(20)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  !Sixth sub-element
                                  ne=(grid_ne-1)*ELEMENT_FACTOR+6
                                  ELEMENT_NODES(1)=np
                                  ELEMENT_NODES(2)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(3)=np+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(4)=np+3+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(5)=np+TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(6)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(7)=np+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(8)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(9)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(10)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(11)=np+TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(12)=np+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(13)=np+1+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(14)=np+2+3*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(15)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(16)=np+2+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(17)=np+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(18)=np+1+TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(19)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+2*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  ELEMENT_NODES(20)=np+1+2*TOTAL_NUMBER_OF_NODES_XI(1)+3*TOTAL_NUMBER_OF_NODES_XI(1)* &
                                    & TOTAL_NUMBER_OF_NODES_XI(2)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(REGULAR_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                      & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                CASE DEFAULT
                                  LOCAL_ERROR="The simplex basis interpolation order of "// &
                                    & TRIM(NUMBER_TO_VSTRING(BASIS%INTERPOLATION_ORDER(1),"*",ERR,ERROR))// &
                                    & " is invalid."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                END SELECT
                              CASE DEFAULT
                                LOCAL_ERROR="The simplex number of xi directions of "// &
                                  & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))// &
                                  & " is invalid."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !ne1
                  ENDDO !ne2
                ENDDO !ne3
                CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,ERR,ERROR,*999)
              ENDDO !basis_idx
              !Finish the mesh
              CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,ERR,ERROR,*999)
            CASE DEFAULT
              CALL FLAG_ERROR("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Bases are not allocated.",ERR,ERROR,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Regular mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)

    CALL EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    CALL ERRORS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start to create the regular generated mesh type
  SUBROUTINE GENERATED_MESH_ELLIPSOID_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH
    TYPE(BASIS_TYPE), POINTER :: BASIS1,BASIS2
    INTEGER(INTG), ALLOCATABLE :: NUMBER_ELEMENTS_XI(:)!,NUMBER_OF_NODES_XIC(:)
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: BASIS_NUMBER_OF_NODES,CORNER_NUMBER_OF_NODES
    INTEGER(INTG) :: ne1,ne2,ne3,nn1,nn2,nn3,from1,from2,from3,nn,ne,mc
    INTEGER(INTG), ALLOCATABLE :: APEX_ELEMENT_NODES(:), WALL_ELEMENT_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: APEX_ELEMENT_NODES_USER_NUMBERS(:), WALL_ELEMENT_NODES_USER_NUMBERS(:)
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),CORNER_NODES(:,:,:),EIDX(:,:,:)
    REAL(DP) :: DELTA(3),DELTAi(3)
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS

    CALL ENTERS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      ELLIPSOID_MESH=>GENERATED_MESH%ELLIPSOID_MESH
      IF(ASSOCIATED(ELLIPSOID_MESH)) THEN
        REGION=>GENERATED_MESH%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
            SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              !Determine the coordinate system and create the regular mesh for that system
              ELLIPSOID_MESH%MESH_DIMENSION=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              NUMBER_OF_DIMENSIONS=ELLIPSOID_MESH%MESH_DIMENSION
              IF(NUMBER_OF_DIMENSIONS==3) THEN ! hard-coded for 3D only
                IF(.NOT.ALLOCATED(ELLIPSOID_MESH%ORIGIN)) THEN
                  ALLOCATE(ELLIPSOID_MESH%ORIGIN(NUMBER_OF_DIMENSIONS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate origin.",ERR,ERROR,*999)
                  ELLIPSOID_MESH%ORIGIN=0.0_DP
                ENDIF
                IF(SIZE(ELLIPSOID_MESH%ORIGIN)==ELLIPSOID_MESH%MESH_DIMENSION) THEN
                  IF(SIZE(ELLIPSOID_MESH%ELLIPSOID_EXTENT)==4) THEN
                    IF(ALLOCATED(ELLIPSOID_MESH%BASES)) THEN
                      IF(MOD(SIZE(ELLIPSOID_MESH%BASES),2)==0) THEN
                        ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                        NUMBER_ELEMENTS_XI=ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI
                        !Calculate total number of nodes from all bases and start mesh
                        CORNER_NUMBER_OF_NODES=NUMBER_ELEMENTS_XI(1)*(NUMBER_ELEMENTS_XI(2)+1)*(NUMBER_ELEMENTS_XI(3)+1)- &
                          & (NUMBER_ELEMENTS_XI(1)-1)*(NUMBER_ELEMENTS_XI(3)+1)
                        TOTAL_NUMBER_OF_NODES=CORNER_NUMBER_OF_NODES
                        DO mc=1,SIZE(ELLIPSOID_MESH%BASES),2
                          BASIS1=>ELLIPSOID_MESH%BASES(mc)%PTR
                          BASIS_NUMBER_OF_NODES=NUMBER_ELEMENTS_XI(1)*(BASIS1%NUMBER_OF_NODES_XIC(1)-1)* &
                            & (NUMBER_ELEMENTS_XI(2)*(BASIS1%NUMBER_OF_NODES_XIC(2)-1)+1)* &
                            & (NUMBER_ELEMENTS_XI(3)*(BASIS1%NUMBER_OF_NODES_XIC(3)-1)+1)- &
                            & (NUMBER_ELEMENTS_XI(1)*(BASIS1%NUMBER_OF_NODES_XIC(1)-1)-1)* &
                            & (NUMBER_ELEMENTS_XI(3)*(BASIS1%NUMBER_OF_NODES_XIC(3)-1)+1)
                          TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES-CORNER_NUMBER_OF_NODES
                        ENDDO
                        NULLIFY(NODES)
                        CALL NODES_CREATE_START(REGION,TOTAL_NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
                        !Finish the nodes creation
                        CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
                        !Create the mesh
                        CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%REGION, &
                          & SIZE(NUMBER_ELEMENTS_XI,1), GENERATED_MESH%MESH,ERR,ERROR,*999)
                        !Create the elements
                        CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
                        DO mc=1,SIZE(ELLIPSOID_MESH%BASES),2
                          IF((ELLIPSOID_MESH%BASES(mc)%PTR%NUMBER_OF_COLLAPSED_XI==0).AND. &
                              & (ELLIPSOID_MESH%BASES(mc+1)%PTR%NUMBER_OF_COLLAPSED_XI>0))THEN
                            !test for collapsed nodes and force non collapsed to wall elements and collapsed to apex elements
                            BASIS1=>ELLIPSOID_MESH%BASES(mc)%PTR
                            BASIS2=>ELLIPSOID_MESH%BASES(mc+1)%PTR
                          ELSE
                            CALL FLAG_ERROR("For each basis, one non collapsed version (basis1) and one collapsed "// &
                                "version (basis2) is needed.",ERR,ERROR,*999)
                          ENDIF
                          SELECT CASE(BASIS1%TYPE)
                          !should also test for basis2
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                            IF(BASIS1%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1).AND. &
                                & BASIS2%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                              IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) &
                                  & CALL FLAG_ERROR("Must have 1 or more elements in all directions.",ERR,ERROR,*999)
                              IF(NUMBER_ELEMENTS_XI(1)<3) &
                                & CALL FLAG_ERROR("Need >2 elements around the circumferential direction.", &
                                & ERR,ERROR,*999)
                              !IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                              !    & CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
                              IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
                              IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
                              IF(ALLOCATED(CORNER_NODES)) DEALLOCATE(CORNER_NODES)
                              CALL GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,BASIS1% &
                                  NUMBER_OF_NODES_XIC, ELLIPSOID_MESH%ELLIPSOID_EXTENT, TOTAL_NUMBER_OF_NODES, &
                                  TOTAL_NUMBER_OF_ELEMENTS, NIDX,CORNER_NODES,EIDX,DELTA,DELTAi,ERR,ERROR,*999)
                              IF(mc==1) THEN
                                CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS, &
                                  & ERR,ERROR,*999)
                              ENDIF

                              !Create the default node set
                              !TODO we finish create after the nodes are initialised?

                              NULLIFY(MESH_ELEMENTS)
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,mc/2+1,BASIS1,MESH_ELEMENTS, &
                                  ERR, ERROR,*999)
                              !Set the elements for the ellipsoid mesh
                              IF(ALLOCATED(WALL_ELEMENT_NODES)) DEALLOCATE(WALL_ELEMENT_NODES)
                              IF(ALLOCATED(APEX_ELEMENT_NODES)) DEALLOCATE(APEX_ELEMENT_NODES)
                              IF(ALLOCATED(WALL_ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(WALL_ELEMENT_NODES_USER_NUMBERS)
                              IF(ALLOCATED(APEX_ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(APEX_ELEMENT_NODES_USER_NUMBERS)
                              ALLOCATE(WALL_ELEMENT_NODES(BASIS1%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate wall element nodes.",ERR,ERROR,*999)
                              ALLOCATE(APEX_ELEMENT_NODES(BASIS2%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate apex element nodes.",ERR,ERROR,*999)
                              ALLOCATE(WALL_ELEMENT_NODES_USER_NUMBERS(BASIS1%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate wall element nodes.",ERR,ERROR,*999)
                              ALLOCATE(APEX_ELEMENT_NODES_USER_NUMBERS(BASIS2%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate apex element nodes.",ERR,ERROR,*999)
                              ! calculate element topology (nodes per each element)
                              ! the idea is to translate given (r,theta,z) to NIDX equivalents, which include interior nodes
                              ne=0
                              nn=0
                              !fromJ=global J direction counting number of first node in element in J direction
                              DO ne3=1,NUMBER_ELEMENTS_XI(3)
                                from3=NINT(DELTA(3)*(ne3-1)/DELTAi(3)+1)
                                ne2=1
                                from2=NINT(DELTA(2)*(ne2-1)/DELTAi(2)+1)
                                !apex elements
                                DO ne1=1,NUMBER_ELEMENTS_XI(1)
                                  from1=NINT(DELTA(1)*(ne1-1)/DELTAi(1)+1)
                                  nn=0
                                  ! number of nodes in an element is dependent on basis used
                                  DO nn3=from3,from3+BASIS2%NUMBER_OF_NODES_XIC(3)-1
                                    nn2=1
                                    nn1=1
                                    !central axis nodes
                                    nn=nn+1
                                    APEX_ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                    DO nn2=from2+1,from2+BASIS2%NUMBER_OF_NODES_XIC(2)-1
                                      DO nn1=from1,from1+BASIS2%NUMBER_OF_NODES_XIC(1)-1
                                        nn=nn+1
                                        ! circumferential loop-around
                                        IF(nn1>SIZE(NIDX,1)) THEN
                                          APEX_ELEMENT_NODES(nn)=NIDX(1,nn2,nn3)
                                        ELSE
                                          APEX_ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                        ENDIF
                                      ENDDO ! nn1
                                    ENDDO ! nn2
                                  ENDDO ! nn3
                                  ne=ne+1
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(ne,MESH_ELEMENTS,BASIS2,ERR,ERROR,*999)
                                  CALL COMPONENT_NODES_TO_USER_NUMBERS(ELLIPSOID_MESH%GENERATED_MESH,mc,APEX_ELEMENT_NODES, &
                                      & APEX_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                    APEX_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                ENDDO ! ne1
                                !wall elements
                                DO ne2=2,NUMBER_ELEMENTS_XI(2)
                                  from2=NINT(DELTA(2)*(ne2-1)/DELTAi(2)+1)
                                  DO ne1=1,NUMBER_ELEMENTS_XI(1)
                                    from1=NINT(DELTA(1)*(ne1-1)/DELTAi(1)+1)
                                    nn=0
                                    ! number of nodes in an element is dependent on basis used
                                    DO nn3=from3,from3+BASIS1%NUMBER_OF_NODES_XIC(3)-1
                                      DO nn2=from2,from2+BASIS1%NUMBER_OF_NODES_XIC(2)-1
                                        DO nn1=from1,from1+BASIS1%NUMBER_OF_NODES_XIC(1)-1
                                          nn=nn+1
                                          ! circumferential loop-around
                                          IF(nn1>SIZE(NIDX,1)) THEN
                                            WALL_ELEMENT_NODES(nn)=NIDX(1,nn2,nn3)
                                          ELSE
                                            WALL_ELEMENT_NODES(nn)=NIDX(nn1,nn2,nn3)
                                          ENDIF
                                        ENDDO ! nn1
                                      ENDDO ! nn2
                                    ENDDO ! nn3
                                    ne=ne+1
                                    CALL COMPONENT_NODES_TO_USER_NUMBERS(ELLIPSOID_MESH%GENERATED_MESH,mc,WALL_ELEMENT_NODES, &
                                        & WALL_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS, &
                                        & WALL_ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                  ENDDO ! ne1
                                ENDDO ! ne2
                              ENDDO ! ne3
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,ERR,ERROR,*999)
                            ELSE
                              CALL FLAG_ERROR("The number of xi directions of the given basis does not match the size of &
                                &the number of elements for the mesh.",ERR,ERROR,*999)
                            ENDIF
                          CASE(BASIS_SIMPLEX_TYPE)
                            CALL FLAG_ERROR("Ellipsoid meshes with simplex basis types is not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            CALL FLAG_ERROR("Basis type is either invalid or not implemented.",ERR,ERROR,*999)
                          END SELECT
                        ENDDO
                        !Finish the mesh
                        CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("An ellipsoid mesh requires a collapsed basis for each basis,"// &
                            & " so there must be n*2 bases.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Bases is not allocated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("For an ellipsoid mesh the following measures need to be given: &
                        & LONG_AXIS, SHORT_AXIS, WALL_THICKNESS and CUTOFF_ANGLE.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                      &the origin.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Ellipsoid mesh requires a 3 dimensional coordinate system.",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Ellipsoid mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    IF(ALLOCATED(CORNER_NODES)) DEALLOCATE(CORNER_NODES)
    IF(ALLOCATED(NUMBER_ELEMENTS_XI)) DEALLOCATE(NUMBER_ELEMENTS_XI)
    IF(ALLOCATED(WALL_ELEMENT_NODES)) DEALLOCATE(WALL_ELEMENT_NODES)
    IF(ALLOCATED(APEX_ELEMENT_NODES)) DEALLOCATE(APEX_ELEMENT_NODES)
    CALL ERRORS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ELLIPSOID_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_CREATE_FINISH
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
    INTEGER(INTG) :: CORNER_NUMBER_OF_NODES,BASIS_NUMBER_OF_NODES
    INTEGER(INTG) :: ne1,ne2,ne3,nn1,nn2,nn3,from1,from2,from3,nn,ne,basis_idx
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:),ELEMENT_NODES_USER_NUMBERS(:)
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
                    IF(ALLOCATED(CYLINDER_MESH%BASES)) THEN
                      ALLOCATE(NUMBER_ELEMENTS_XI(SIZE(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of elements xi.",ERR,ERROR,*999)
                      NUMBER_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
                      CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%REGION,SIZE(NUMBER_ELEMENTS_XI,1), &
                        & GENERATED_MESH%MESH,ERR,ERROR,*999)
                      CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(CYLINDER_MESH%BASES),ERR,ERROR,*999)
                      !Calculate number of nodes
                      CORNER_NUMBER_OF_NODES=(NUMBER_ELEMENTS_XI(3)+1)*NUMBER_ELEMENTS_XI(2)*(NUMBER_ELEMENTS_XI(1)+1)
                      TOTAL_NUMBER_OF_NODES=CORNER_NUMBER_OF_NODES
                      DO basis_idx=1,SIZE(CYLINDER_MESH%BASES)
                        BASIS=>CYLINDER_MESH%BASES(basis_idx)%PTR
                        IF(ASSOCIATED(BASIS)) THEN
                          BASIS_NUMBER_OF_NODES=((BASIS%NUMBER_OF_NODES_XIC(3)-1)*NUMBER_ELEMENTS_XI(3)+1)* &
                              & ((BASIS%NUMBER_OF_NODES_XIC(2)-1)*NUMBER_ELEMENTS_XI(2))* &
                              & ((BASIS%NUMBER_OF_NODES_XIC(1)-1)*NUMBER_ELEMENTS_XI(1)+1)
                          TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+BASIS_NUMBER_OF_NODES-CORNER_NUMBER_OF_NODES
                        ELSE
                          CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO
                      NULLIFY(NODES)
                      CALL NODES_CREATE_START(REGION,TOTAL_NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
                      !Finish the nodes creation
                      CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
                      !Set the total number of elements
                      TOTAL_NUMBER_OF_ELEMENTS=NUMBER_ELEMENTS_XI(1)*NUMBER_ELEMENTS_XI(2)*NUMBER_ELEMENTS_XI(3)
                      CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
                      DO basis_idx=1,SIZE(CYLINDER_MESH%BASES)
                        BASIS=>CYLINDER_MESH%BASES(basis_idx)%PTR
                        IF(ASSOCIATED(BASIS)) THEN
                          SELECT CASE(BASIS%TYPE)
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                            IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                              IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) &
                                & CALL FLAG_ERROR("Must have 1 or more elements in all directions.",ERR,ERROR,*999)
                              IF(NUMBER_ELEMENTS_XI(2)<3) &
                                CALL FLAG_ERROR("Need >2 elements around the circumferential direction.",ERR,ERROR,*999)
                              IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED))  &
                                & CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented.",ERR,ERROR,*999)
                              !Calculate nodes and element sizes
                              IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
                              IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
                              CALL GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,BASIS%NUMBER_OF_NODES_XIC, &
                                & CYLINDER_MESH%CYLINDER_EXTENT, TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS, &
                                & NIDX,EIDX,DELTA,DELTAi,ERR,ERROR,*999)
                              !Set the elements for the cylinder mesh
                              IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
                              IF(ALLOCATED(ELEMENT_NODES_USER_NUMBERS)) DEALLOCATE(ELEMENT_NODES_USER_NUMBERS)
                              ALLOCATE(ELEMENT_NODES_USER_NUMBERS(BASIS%NUMBER_OF_NODES),STAT=ERR)
                              ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element nodes.",ERR,ERROR,*999)
                              !Create the elements
                              NULLIFY(MESH_ELEMENTS)
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(GENERATED_MESH%MESH,basis_idx,BASIS,MESH_ELEMENTS, &
                                  & ERR,ERROR,*999)
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
                                    DO nn3=from3,from3+BASIS%NUMBER_OF_NODES_XIC(3)-1
                                      DO nn2=from2,from2+BASIS%NUMBER_OF_NODES_XIC(2)-1
                                        DO nn1=from1,from1+BASIS%NUMBER_OF_NODES_XIC(1)-1
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
                                    CALL COMPONENT_NODES_TO_USER_NUMBERS(CYLINDER_MESH%GENERATED_MESH,basis_idx,ELEMENT_NODES, &
                                        & ELEMENT_NODES_USER_NUMBERS,ERR,ERROR,*999)
                                    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS,ELEMENT_NODES_USER_NUMBERS, &
                                        & ERR,ERROR,*999)
                                  ENDDO ! ne1
                                ENDDO ! ne2
                              ENDDO ! ne3
                              CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH_ELEMENTS,ERR,ERROR,*999)
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
                      ENDDO
                      !Finish the mesh
                      CALL MESH_CREATE_FINISH(GENERATED_MESH%MESH,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Bases are not allocated.",ERR,ERROR,*999)
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
      IF(ALLOCATED(CYLINDER_MESH%BASES)) DEALLOCATE(CYLINDER_MESH%BASES)
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

  !>Initialise the cylinder generated mesh type
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
      IF(ALLOCATED(REGULAR_MESH%BASE_VECTORS)) DEALLOCATE(REGULAR_MESH%BASE_VECTORS)
      IF(ALLOCATED(REGULAR_MESH%BASES)) DEALLOCATE(REGULAR_MESH%BASES)
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

  !>Finalise ellipsoid mesh type
  SUBROUTINE GENERATED_MESH_ELLIPSOID_FINALISE(ELLIPSOID_MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_ELLIPSOID_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELLIPSOID_MESH)) THEN
      IF(ALLOCATED(ELLIPSOID_MESH%ORIGIN)) DEALLOCATE(ELLIPSOID_MESH%ORIGIN)
      IF(ALLOCATED(ELLIPSOID_MESH%ELLIPSOID_EXTENT)) DEALLOCATE(ELLIPSOID_MESH%ELLIPSOID_EXTENT)
      IF(ALLOCATED(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI)) DEALLOCATE(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI)
      IF(ALLOCATED(ELLIPSOID_MESH%BASES)) DEALLOCATE(ELLIPSOID_MESH%BASES)
      DEALLOCATE(ELLIPSOID_MESH)
    ENDIF

    CALL EXITS("GENERATED_MESH_ELLIPSOID_FINALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL ERRORS("GENERATED_MESH_ELLIPSOID_FINALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ELLIPSOID_FINALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the ellipsoid generated mesh type
  SUBROUTINE GENERATED_MESH_ELLIPSOID_INITIALISE(GENERATED_MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("GENERATED_MESH_ELLIPSOID_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
        CALL FLAG_ERROR("Ellipsoid mesh is already associated for this generated mesh.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GENERATED_MESH%ELLIPSOID_MESH,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ellipsoid generated mesh.",ERR,ERROR,*999)
        GENERATED_MESH%ELLIPSOID_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_ELLIPSOID_MESH_TYPE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("GENERATED_MESH_ELLIPSOID_INITIALISE")
    RETURN
999 CALL GENERATED_MESH_ELLIPSOID_FINALISE(GENERATED_MESH%ELLIPSOID_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GENERATED_MESH_ELLIPSOID_INITIALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ELLIPSOID_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_INITIALISE

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
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
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
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_FINALISE(GENERATED_MESH%ELLIPSOID_MESH,ERR,ERROR,*999)
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
            CALL GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE(GENERATED_MESH%REGULAR_MESH, &
                & FIELD,ERR,ERROR,*999)
          CASE(GENERATED_MESH_POLAR_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
            CALL GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE(GENERATED_MESH%CYLINDER_MESH, &
                & FIELD,ERR,ERROR,*999)
          CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
            CALL GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE(GENERATED_MESH%ELLIPSOID_MESH, &
                & FIELD,ERR,ERROR,*999)
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

  !>Returns the region for a generated mesh accounting for regions and interfaces
  SUBROUTINE GENERATED_MESH_REGION_GET(GENERATED_MESH,REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the region for
    TYPE(REGION_TYPE), POINTER :: REGION !<On return, the generated meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_REGION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(REGION)) THEN
        CALL FLAG_ERROR("Region is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(REGION)
        NULLIFY(INTERFACE)
        REGION=>GENERATED_MESH%REGION
        IF(.NOT.ASSOCIATED(REGION)) THEN
          INTERFACE=>GENERATED_MESH%INTERFACE
          IF(ASSOCIATED(INTERFACE)) THEN
            PARENT_REGION=>INTERFACE%PARENT_REGION
            IF(ASSOCIATED(PARENT_REGION)) THEN
              REGION=>PARENT_REGION
            ELSE
              LOCAL_ERROR="The parent region not associated for generated mesh number "// &
                & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//" of interface number "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The region or interface is not associated for generated mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_REGION_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_REGION_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGION_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGION_GET

  !
  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the regular mesh. Any derivative values for the nodes are calculated from an average straight line approximation.
  SUBROUTINE GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE(REGULAR_MESH,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH !<A pointer to the regular mesh object
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local variables
    INTEGER(INTG) :: component_idx,component_idx2,derivative_idx,derivative1,derivative2, &
      & DERIVATIVES_NUMBER_OF_LINES(MAXIMUM_GLOBAL_DERIV_NUMBER),dof,global_node,global_node1,global_node2, &
      & component_node,component_node1,component_node2,MESH_COMPONENT, &
      & line,local_line_idx,node_idx,node1,node2,node_position_idx(3),node_position_idx2(3), &
      & partial_derivative,TOTAL_NUMBER_OF_NODES_XI(3),xi_idx
    REAL(DP) :: DELTA(MAXIMUM_GLOBAL_DERIV_NUMBER),DELTA_COORD(3,3),LENGTH,MY_ORIGIN(3),MY_EXTENT(3),VALUE,VECTOR(3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    LOGICAL :: HAVE_DERIVATIVES
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_REGULAR_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGULAR_MESH)) THEN
      IF(ASSOCIATED(FIELD)) THEN
        NULLIFY(COORDINATE_SYSTEM)
        CALL FIELD_COORDINATE_SYSTEM_GET(FIELD,COORDINATE_SYSTEM,ERR,ERROR,*999)
        IF(COORDINATE_SYSTEM%TYPE==COORDINATE_RECTANGULAR_CARTESIAN_TYPE) THEN

          MY_ORIGIN=0.0_DP
          MY_EXTENT=0.0_DP
          MY_ORIGIN(1:REGULAR_MESH%COORDINATE_DIMENSION)=REGULAR_MESH%ORIGIN(1:REGULAR_MESH%COORDINATE_DIMENSION)
          MY_EXTENT(1:REGULAR_MESH%COORDINATE_DIMENSION)=REGULAR_MESH%MAXIMUM_EXTENT(1:REGULAR_MESH%COORDINATE_DIMENSION)
          DELTA_COORD=0.0_DP
          TOTAL_NUMBER_OF_NODES_XI=1
          IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              HAVE_DERIVATIVES=.FALSE.
              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                FIELD_VARIABLE_COMPONENT=>FIELD_VARIABLE%COMPONENTS(component_idx)
                MESH_COMPONENT=FIELD_VARIABLE_COMPONENT%MESH_COMPONENT_NUMBER
                IF(FIELD_VARIABLE_COMPONENT%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                  DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                    TOTAL_NUMBER_OF_NODES_XI(xi_idx)=(REGULAR_MESH%BASES(MESH_COMPONENT)% &
                        & PTR%NUMBER_OF_NODES_XIC(xi_idx)-2)*REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(xi_idx)+ &
                        & REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(xi_idx)+1
                  ENDDO !xi_idx
                  DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                    DELTA_COORD(1:REGULAR_MESH%COORDINATE_DIMENSION,xi_idx)= &
                      & REGULAR_MESH%BASE_VECTORS(1:REGULAR_MESH%COORDINATE_DIMENSION,xi_idx)/ &
                      & REAL(TOTAL_NUMBER_OF_NODES_XI(xi_idx)-1,DP)
                  ENDDO !xi_idx
                  !Update geometric parameters in this computational domain only
                  DOMAIN=>FIELD_VARIABLE_COMPONENT%DOMAIN
                  DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                  DOMAIN_LINES=>DOMAIN%TOPOLOGY%LINES
                  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                    !Domain nodes are only the nodes for this mesh component
                    global_node=DOMAIN_NODES%NODES(node_idx)%GLOBAL_NUMBER
                    component_node=USER_NUMBER_TO_COMPONENT_NODE(REGULAR_MESH%GENERATED_MESH, &
                        & MESH_COMPONENT,global_node,ERR,ERROR)
                    node_position_idx(3)=(component_node-1)/(TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))+1
                    node_position_idx(2)=MOD(component_node-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))/ &
                      & TOTAL_NUMBER_OF_NODES_XI(1)+1
                    node_position_idx(1)=MOD(MOD(component_node-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1)), &
                      & TOTAL_NUMBER_OF_NODES_XI(1))+1
                    dof=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                    VALUE=0.0_DP
                    DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                      VALUE=VALUE+REAL(node_position_idx(xi_idx)-1,DP)*DELTA_COORD(component_idx,xi_idx)
                    ENDDO !xi_idx
                    VALUE=MY_ORIGIN(component_idx)+VALUE
                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof,VALUE, &
                      & ERR,ERROR,*999)
                    !Calculate derivatives
                    IF(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES>1) THEN
                      HAVE_DERIVATIVES=.TRUE.
                      DERIVATIVES_NUMBER_OF_LINES=0
                      DELTA=0.0_DP
                      DO local_line_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES
                        line=DOMAIN_NODES%NODES(node_idx)%NODE_LINES(local_line_idx)
                        node1=DOMAIN_LINES%LINES(line)%NODES_IN_LINE(1)
                        global_node1=DOMAIN_NODES%NODES(node1)%GLOBAL_NUMBER
                        component_node1=USER_NUMBER_TO_COMPONENT_NODE(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                            & global_node1,ERR,ERROR)
                        node2=DOMAIN_LINES%LINES(line)%NODES_IN_LINE(DOMAIN_LINES%LINES(line)%BASIS%NUMBER_OF_NODES)
                        global_node2=DOMAIN_NODES%NODES(node2)%GLOBAL_NUMBER
                        component_node2=USER_NUMBER_TO_COMPONENT_NODE(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                            & global_node2,ERR,ERROR)
                        derivative1=DOMAIN_LINES%LINES(line)%DERIVATIVES_IN_LINE(2,1)
                        derivative2=DOMAIN_LINES%LINES(line)%DERIVATIVES_IN_LINE(2,DOMAIN_LINES%LINES(line)%BASIS%NUMBER_OF_NODES)
                        VALUE=0.0_DP
                        IF(node1==node_idx) THEN
                          node_position_idx2(3)=(component_node2-1)/(TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))+1
                          node_position_idx2(2)=MOD(component_node2-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))/ &
                            & TOTAL_NUMBER_OF_NODES_XI(1)+1
                          node_position_idx2(1)=MOD(MOD(component_node2-1,TOTAL_NUMBER_OF_NODES_XI(2)* &
                            & TOTAL_NUMBER_OF_NODES_XI(1)),TOTAL_NUMBER_OF_NODES_XI(1))+1
                          DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                            VALUE=VALUE+REAL(node_position_idx2(xi_idx)-node_position_idx(xi_idx),DP)* &
                              & DELTA_COORD(component_idx,xi_idx)
                          ENDDO !xi_idx
                          DERIVATIVES_NUMBER_OF_LINES(derivative1)=DERIVATIVES_NUMBER_OF_LINES(derivative1)+1
                          DELTA(derivative1)=DELTA(derivative1)+VALUE
                        ELSE
                          node_position_idx2(3)=(component_node1-1)/(TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))+1
                          node_position_idx2(2)=MOD(component_node1-1,TOTAL_NUMBER_OF_NODES_XI(2)*TOTAL_NUMBER_OF_NODES_XI(1))/ &
                            & TOTAL_NUMBER_OF_NODES_XI(1)+1
                          node_position_idx2(1)=MOD(MOD(component_node1-1,TOTAL_NUMBER_OF_NODES_XI(2)* &
                            & TOTAL_NUMBER_OF_NODES_XI(1)),TOTAL_NUMBER_OF_NODES_XI(1))+1
                          DO xi_idx=1,REGULAR_MESH%MESH_DIMENSION
                            VALUE=VALUE+REAL(node_position_idx(xi_idx)-node_position_idx2(xi_idx),DP)* &
                              & DELTA_COORD(component_idx,xi_idx)
                          ENDDO !xi_idx
                          DERIVATIVES_NUMBER_OF_LINES(derivative2)=DERIVATIVES_NUMBER_OF_LINES(derivative2)+1
                          DELTA(derivative2)=DELTA(derivative2)+VALUE
                        ENDIF
                      ENDDO !local_line_idx
                      DO derivative_idx=1,MAXIMUM_GLOBAL_DERIV_NUMBER
                        IF(DERIVATIVES_NUMBER_OF_LINES(derivative_idx)>0) THEN
                          DELTA(derivative_idx)=DELTA(derivative_idx)/REAL(DERIVATIVES_NUMBER_OF_LINES(derivative_idx),DP)
                          dof=FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative_idx,node_idx)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & dof,DELTA(derivative_idx),ERR,ERROR,*999)
                        ENDIF
                      ENDDO !derivative_idx
                    ENDIF
                  ENDDO !node_idx
                ELSE
                  LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                    & " does not have node based interpolation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
              IF(HAVE_DERIVATIVES) THEN
                !Normalise the arclength derivative vectors.
                NULLIFY(GEOMETRIC_PARAMETERS)
                CALL FIELD_PARAMETER_SET_DATA_GET(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
                  & ERR,ERROR,*999)
!!TODO: Don't loop over all components somehow????
                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                  DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                  DOMAIN_LINES=>DOMAIN%TOPOLOGY%LINES
                  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                    DO derivative_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                      partial_derivative=DOMAIN_NODES%NODES(node_idx)%PARTIAL_DERIVATIVE_INDEX(derivative_idx)
                      IF(partial_derivative==PART_DERIV_S1.OR.partial_derivative==PART_DERIV_S2.OR. &
                        & partial_derivative==PART_DERIV_S3) THEN
                        LENGTH=0.0_DP
                        VECTOR=0.0_DP
                        DO component_idx2=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                          dof=FIELD_VARIABLE%COMPONENTS(component_idx2)%PARAM_TO_DOF_MAP% &
                            & NODE_PARAM2DOF_MAP(derivative_idx,node_idx)
                          VECTOR(component_idx2)=GEOMETRIC_PARAMETERS(dof)
                          LENGTH=LENGTH+VECTOR(component_idx2)**2
                        ENDDO !component_idx2
                        LENGTH=SQRT(LENGTH)
                        IF(LENGTH>ZERO_TOLERANCE) THEN
                          VECTOR=VECTOR/LENGTH
                          DO component_idx2=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            dof=FIELD_VARIABLE%COMPONENTS(component_idx2)%PARAM_TO_DOF_MAP% &
                              & NODE_PARAM2DOF_MAP(derivative_idx,node_idx)
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof, &
                              & VECTOR(component_idx2),ERR,ERROR,*999)
                          ENDDO !component_idx2
                        ENDIF
                      ENDIF
                    ENDDO !derivative_idx
                  ENDDO !node_idx
                ENDDO !component_idx
                CALL FIELD_PARAMETER_SET_DATA_RESTORE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
                  & ERR,ERROR,*999)
              ENDIF
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
        ELSE
          CALL FLAG_ERROR("Non rectangular Cartesian coordinate systems are not implemented.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Regular mesh is not associated.",ERR,ERROR,*999)
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
    INTEGER(INTG) :: NUMBER_ELEMENTS_XI(3),NUMBER_OF_NODES_XIC(3)
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3),INTERPOLATION_TYPES(3)
    INTEGER(INTG) :: component_idx,xi_idx
    INTEGER(INTG) :: np,global_np,component_np,ny,nk
    INTEGER(INTG) :: NUMBER_OF_PLANAR_NODES,SCALING_TYPE,MESH_COMPONENT
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: node_idx(3) ! holds r,theta,z indices
    REAL(DP) :: DELTA(3),DELTAi(3),POLAR_COORDS(3),RECT_COORDS(3)
    REAL(DP) :: CYLINDER_EXTENT(3),DERIV
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_CYLINDER_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    ! assign to the field
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
          CALL FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,ERR,ERROR,*999)
          IF(SCALING_TYPE/=FIELD_UNIT_SCALING) &
            & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"  Note: If the cylinder looks wonky, set field scaling to&
            & unit scaling type.",ERR,ERROR,*999)
          DO component_idx=1,3
            INTERPOLATION_TYPES(component_idx)=FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE
          ENDDO
          IF(ALL(INTERPOLATION_TYPES==FIELD_NODE_BASED_INTERPOLATION)) THEN
            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              FIELD_VARIABLE_COMPONENT=>FIELD_VARIABLE%COMPONENTS(component_idx)
              MESH_COMPONENT=FIELD_VARIABLE_COMPONENT%MESH_COMPONENT_NUMBER
              ! calculate the total number of nodes in each xi direction
              BASIS=>CYLINDER_MESH%BASES(MESH_COMPONENT)%PTR
              NUMBER_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
              NUMBER_OF_NODES_XIC=BASIS%NUMBER_OF_NODES_XIC
              DO xi_idx=1,3
                TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
              ENDDO
              TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
              NUMBER_OF_PLANAR_NODES=TOTAL_NUMBER_NODES_XI(1)*TOTAL_NUMBER_NODES_XI(2)
              DOMAIN=>FIELD_VARIABLE%COMPONENTS(MESH_COMPONENT)%DOMAIN
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              ! calculate DELTAi now
              CYLINDER_EXTENT=CYLINDER_MESH%CYLINDER_EXTENT
              DELTA(1)=(CYLINDER_EXTENT(2)-CYLINDER_EXTENT(1))/NUMBER_ELEMENTS_XI(1)
              DELTA(2)=TWOPI/NUMBER_ELEMENTS_XI(2)
              DELTA(3)=CYLINDER_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
              DO xi_idx=1,3
                DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
              ENDDO
              DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
                global_np=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
                component_np=USER_NUMBER_TO_COMPONENT_NODE(CYLINDER_MESH%GENERATED_MESH, &
                    & MESH_COMPONENT,global_np,ERR,ERROR)
                ! calculate node_idx which will be used to calculate (r,theta,z) then (x,y,z)
                component_np=component_np-1 ! let's go 0-based index for a bit
                node_idx(3)=component_np/NUMBER_OF_PLANAR_NODES
                node_idx(2)=(component_np-(node_idx(3))*NUMBER_OF_PLANAR_NODES)/TOTAL_NUMBER_NODES_XI(1)
                node_idx(1)=MOD(component_np-(node_idx(3))*NUMBER_OF_PLANAR_NODES,TOTAL_NUMBER_NODES_XI(1))
                DO xi_idx=1,3
                  POLAR_COORDS(xi_idx)=node_idx(xi_idx)*DELTAi(xi_idx)
                ENDDO
                POLAR_COORDS(1)=node_idx(1)*DELTAi(1)+CYLINDER_EXTENT(1) ! add the inner radius
                RECT_COORDS(1)=POLAR_COORDS(1)*COS(POLAR_COORDS(2))
                RECT_COORDS(2)=POLAR_COORDS(1)*SIN(POLAR_COORDS(2))
                RECT_COORDS(3)=POLAR_COORDS(3)
                RECT_COORDS=RECT_COORDS+CYLINDER_MESH%ORIGIN
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
              ENDDO !np
            ENDDO !component_idx
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

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE(ELLIPSOID_MESH,FIELD,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH !<A pointer to the ellipsoid mesh object
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_TYPE),POINTER :: DOMAIN
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_VARIABLE_COMPONENT
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE,DOMAIN_NUMBER,MESH_COMPONENT,basis_idx
    INTEGER(INTG) :: NUMBER_ELEMENTS_XI(3),NUMBER_OF_NODES_XIC(3)
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3),INTERPOLATION_TYPES(3)
    INTEGER(INTG) :: component_idx,xi_idx
    INTEGER(INTG) :: np,npg,i,j,k, local_node
    INTEGER(INTG) :: SCALING_TYPE!,NUMBER_OF_PLANAR_NODES
    INTEGER(INTG), ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    !INTEGER(INTG) :: node_idx(3) ! holds r,theta,z indices
    REAL(DP) :: DELTA(3),DELTAi(3),RECT_COORDS(3),t,phi,alpha,xi,nu,x,y,z
    REAL(DP) :: ELLIPSOID_EXTENT(4)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(BASIS,DOMAIN,DECOMPOSITION,DOMAIN_NODES,FIELD_VARIABLE,FIELD_VARIABLE_COMPONENT)

    CALL ENTERS("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)

    ! assign to the field
    np=0
    IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
       FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
       IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==3) THEN
             MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
             DO component_idx=2,3
                IF(FIELD_VARIABLE%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER/=MESH_COMPONENT) THEN
                   CALL FLAG_ERROR("Multiple mesh components for geometric components is not implemented.",ERR,ERROR,*999)
                ENDIF
             ENDDO
             basis_idx=MESH_COMPONENT*2-1
             !< Ellipsoid_extent= inner long axis, inner short axis, wall thickness, top angle (from 0)
             ! calculate the total number of nodes in each xi direction
             IF(ALLOCATED(ELLIPSOID_MESH%BASES)) THEN
                !Check that the all geometric bases use the same mesh component
                BASIS=>ELLIPSOID_MESH%BASES(basis_idx)%PTR
                NUMBER_ELEMENTS_XI=ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI
                NUMBER_OF_NODES_XIC=BASIS%NUMBER_OF_NODES_XIC
                DO xi_idx=1,3
                   TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
                ENDDO
                TOTAL_NUMBER_NODES_XI(1)=TOTAL_NUMBER_NODES_XI(1)-1 ! theta loops around so slightly different
                ! calculate DELTAi now
                ELLIPSOID_EXTENT=ELLIPSOID_MESH%ELLIPSOID_EXTENT
                DELTA(1)=TWOPI/NUMBER_ELEMENTS_XI(1)
                DELTA(2)=(PI-ELLIPSOID_EXTENT(4))/NUMBER_ELEMENTS_XI(2)
                DELTA(3)=ELLIPSOID_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
                DO xi_idx=1,3
                   DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
                ENDDO
             ELSE
                CALL FLAG_ERROR("Ellipsoid mesh does not have bases allocated.",ERR,ERROR,*999)
             ENDIF
             CALL FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,ERR,ERROR,*999)
             IF(SCALING_TYPE/=FIELD_UNIT_SCALING) &
                  & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"  Note: If the ellipsoid looks wonky, set field scaling to &
                  & unit scaling type.",ERR,ERROR,*999)
             ! NUMBER_OF_PLANAR_NODES=TOTAL_NUMBER_NODES_XI(1)*TOTAL_NUMBER_NODES_XI(2)
             DO component_idx=1,3
                INTERPOLATION_TYPES(component_idx)=FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE
             ENDDO
             IF(ALL(INTERPOLATION_TYPES==FIELD_NODE_BASED_INTERPOLATION)) THEN
                DOMAIN=>FIELD_VARIABLE%COMPONENTS(1)%DOMAIN ! just grab the first one
                DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                !DECOMPOSITION=>DOMAIN%DECOMPOSITION !\todo: test all these pointers
                DECOMPOSITION=>FIELD%DECOMPOSITION !\todo: test all these pointers
                IF (ELLIPSOID_EXTENT(1)>ELLIPSOID_EXTENT(2)) THEN
                   !Prolate spheroid
                   k=1
                   !inner surface
                   alpha=sqrt((ELLIPSOID_EXTENT(1))**2-(ELLIPSOID_EXTENT(2))**2)
                   !xi=log(ELLIPSOID_EXTENT(1)/alpha+sqrt((ELLIPSOID_EXTENT(1)/alpha)**2+1))
                   xi=acosh(ELLIPSOID_EXTENT(1)/alpha)

                   j=1
                   !apex node
                   np=1
                   npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                   CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                   IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                         CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                              & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                         local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                         IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                         ENDIF !derivatives
                      ENDDO
                   ENDIF

                   DO j=2,TOTAL_NUMBER_NODES_XI(2)
                      !longitudinal loop
                      nu=PI-DELTAi(2)*(j-1)
                      DO i=1,TOTAL_NUMBER_NODES_XI(1)
                         !circumferential loop
                         phi=DELTAi(1)*(i-1)
                         RECT_COORDS(1)=alpha*(sinh(xi)*sin(nu)*cos(phi))
                         RECT_COORDS(2)=alpha*(sinh(xi)*sin(nu)*sin(phi))
                         RECT_COORDS(3)=alpha*(cosh(xi)*cos(nu))
                         np=np+1
                         npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                         CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                         IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                               CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                    & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                               local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                               IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                  CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                               ENDIF !derivatives
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO

                   DO k=2,TOTAL_NUMBER_NODES_XI(3)
                      !transmural loop
                      j=1
                      !apex nodes
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)-(k-1)*(DELTAi(3))
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                            ENDIF !derivatives
                         ENDDO
                      ENDIF

                      DO j=2,TOTAL_NUMBER_NODES_XI(2)
                         !longitudinal loop
                         nu=PI-DELTAi(2)*(j-1)
                         DO i=1,TOTAL_NUMBER_NODES_XI(1)
                            !circumferential loop
                            phi=DELTAi(1)*(i-1)
                            x=alpha*(sinh(xi)*sin(nu)*cos(phi))
                            y=alpha*(sinh(xi)*sin(nu)*sin(phi))
                            z=alpha*(cosh(xi)*cos(nu))
                            !Normal vector from inner surface with length DELTAi(3)(k-1)
                            ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
                            !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
                            t=(DELTAi(3)*(k-1))/sqrt((4*x**2/(ELLIPSOID_EXTENT(2))**4)+ &
                                 & (4*y**2/(ELLIPSOID_EXTENT(2))**4)+(4*z**2/(ELLIPSOID_EXTENT(1))**4))
                            RECT_COORDS(1)=x*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(2)=y*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(3)=z*(1+2*t/(ELLIPSOID_EXTENT(1))**2)
                            np=np+1
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSEIF (ELLIPSOID_EXTENT(1)==ELLIPSOID_EXTENT(2)) THEN
                   !Sphere
                   np=0
                   DO k=1,TOTAL_NUMBER_NODES_XI(3)
                      !transmural loop
                      alpha=ELLIPSOID_EXTENT(1)+(k-1)*(DELTAi(3))
                      j=1
                      !apex nodes
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-alpha
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                            ENDIF !derivatives
                         ENDDO
                      ENDIF

                      DO j=2,TOTAL_NUMBER_NODES_XI(2)
                         !longitudinal loop
                         nu=PI-DELTAi(2)*(j-1)
                         DO i=1,TOTAL_NUMBER_NODES_XI(1)
                            !circumferential loop
                            phi=DELTAi(1)*(i-1)
                            RECT_COORDS(1)=alpha*sin(nu)*cos(phi)
                            RECT_COORDS(2)=alpha*sin(nu)*sin(phi)
                            RECT_COORDS(3)=alpha*cos(nu)
                            np=np+1
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

                ELSEIF (ELLIPSOID_EXTENT(1)<ELLIPSOID_EXTENT(2)) THEN 
                   !Oblate spheroid
                   k=1
                   !inner surface
                   alpha=sqrt((ELLIPSOID_EXTENT(2))**2-(ELLIPSOID_EXTENT(1))**2)
                   !xi=log(ELLIPSOID_EXTENT(1)/alpha+sqrt((ELLIPSOID_EXTENT(1)/alpha)**2+1))
                   xi=acosh(ELLIPSOID_EXTENT(2)/alpha)

                   j=1
                   !apex node
                   np=1
                   npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                   CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                   IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                         CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                              & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                         local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                         IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                            CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                         ENDIF !derivatives
                      ENDDO
                   ENDIF

                   DO j=2,TOTAL_NUMBER_NODES_XI(2)
                      !longitudinal loop
                      nu=-PI/2+DELTAi(2)*(j-1)
                      DO i=1,TOTAL_NUMBER_NODES_XI(1)
                         !circumferential loop
                         phi=DELTAi(1)*(i-1)
                         RECT_COORDS(1)=alpha*(cosh(xi)*cos(nu)*cos(phi))
                         RECT_COORDS(2)=alpha*(cosh(xi)*cos(nu)*sin(phi))
                         RECT_COORDS(3)=alpha*(sinh(xi)*sin(nu))
                         np=np+1
                         npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                         CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                         IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                               CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                    & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                               local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                               IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                  CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                               ENDIF !derivatives
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO

                   DO k=2,TOTAL_NUMBER_NODES_XI(3)
                      !transmural loop
                      j=1
                      !apex nodes
                      RECT_COORDS(1)=0
                      RECT_COORDS(2)=0
                      RECT_COORDS(3)=-ELLIPSOID_EXTENT(1)-(k-1)*(DELTAi(3))
                      np=np+1
                      npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                         DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                 & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                            local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                            IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                               CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                            ENDIF !derivatives
                         ENDDO
                      ENDIF

                      DO j=2,TOTAL_NUMBER_NODES_XI(2)
                         !longitudinal loop
                         nu=-PI/2+DELTAi(2)*(j-1)
                         DO i=1,TOTAL_NUMBER_NODES_XI(1)
                            !circumferential loop
                            phi=DELTAi(1)*(i-1)
                            x=alpha*(cosh(xi)*cos(nu)*cos(phi))
                            y=alpha*(cosh(xi)*cos(nu)*sin(phi))
                            z=alpha*(sinh(xi)*sin(nu))
                            !Normal vector from inner surface with length DELTAi(3)(k-1)
                            ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
                            !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
                            t=(DELTAi(3)*(k-1))/sqrt((4*x**2/(ELLIPSOID_EXTENT(2))**4)+ &
                                 & (4*y**2/(ELLIPSOID_EXTENT(2))**4)+(4*z**2/(ELLIPSOID_EXTENT(1))**4))
                            RECT_COORDS(1)=x*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(2)=y*(1+2*t/(ELLIPSOID_EXTENT(2))**2)
                            RECT_COORDS(3)=z*(1+2*t/(ELLIPSOID_EXTENT(1))**2)
                            np=np+1
                            npg=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,basis_idx,np,ERR,ERROR)
                            CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,npg,MESH_COMPONENT,DOMAIN_NUMBER,ERR,ERROR,*999)
                            IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE) THEN
                               DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,npg, &
                                       & component_idx,RECT_COORDS(component_idx),ERR,ERROR,*999)
                                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(npg)%local_number(1)
                                  IF(DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES>1) THEN
                                     CALL FLAG_ERROR("Not generalized to hermittean elements.",ERR,ERROR,*999)
                                  ENDIF !derivatives
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE 
                   CALL FLAG_ERROR("Not valid long axis - short axis relation",ERR,ERROR,*999)
                ENDIF
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

    CALL EXITS("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN
999 IF(ALLOCATED(NIDX)) DEALLOCATE(NIDX)
    IF(ALLOCATED(EIDX)) DEALLOCATE(EIDX)
    CALL ERRORS("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN 1

  END SUBROUTINE GENERATED_MESH_ELLIPSOID_GEOMETRIC_PARAMETERS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_REGULAR_SURFACE_GET(REGULAR_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH !<A pointer to the regular mesh object
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GENERATED_MESH_ROUTINES_GeneratedMeshRegularSurfaces,GENERATED_MESH_ROUTINES
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: num_dims,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NODE_USER_NUMBER
    REAL(DP) :: DELTA(3),DELTAI(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: node_counter,i,j,k

    CALL ENTERS("GENERATED_MESH_REGULAR_SURFACE_GET",ERR,ERROR,*999)

    IF(ALLOCATED(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
      num_dims=SIZE(REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)
      IF(num_dims==2) THEN
        NUMBER_OF_ELEMENTS_XI(1:2)=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1:2)
        NUMBER_OF_ELEMENTS_XI(3)=1
      ELSE
        NUMBER_OF_ELEMENTS_XI=REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
      ENDIF
      IF(ASSOCIATED(REGULAR_MESH%BASES(MESH_COMPONENT)%PTR)) THEN
        BASIS=>REGULAR_MESH%BASES(MESH_COMPONENT)%PTR
        IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
          !Node that simplex bases have an extra area coordinate so size of number_of_nodes_xic=num_dims+1
          NUMBER_OF_NODES_XI(1:num_dims)=BASIS%NUMBER_OF_NODES_XIC(1:num_dims)
          NUMBER_OF_NODES_XI(num_dims+1:3)=1

          ! build indices first (some of these are dummy arguments)
          CALL GENERATED_MESH_REGULAR_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
              & REGULAR_MESH%MAXIMUM_EXTENT,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAI,ERR,ERROR,*999)
          node_counter=0
          SELECT CASE(SURFACE_TYPE)
          CASE(GENERATED_MESH_REGULAR_LEFT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(1,j,k),ERR,ERROR)
                SURFACE_NODES(node_counter)=NODE_USER_NUMBER
              ENDDO
            ENDDO
            NORMAL_XI=-1
          CASE(GENERATED_MESH_REGULAR_RIGHT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(SIZE(NIDX,1),j,k),ERR,ERROR)
                SURFACE_NODES(node_counter)=NODE_USER_NUMBER
              ENDDO
            ENDDO
            NORMAL_XI=1
          CASE(GENERATED_MESH_REGULAR_TOP_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
                SURFACE_NODES(node_counter)=NODE_USER_NUMBER
              ENDDO
            ENDDO
            NORMAL_XI=3
          CASE(GENERATED_MESH_REGULAR_BOTTOM_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,1),ERR,ERROR)
                SURFACE_NODES(node_counter)=NODE_USER_NUMBER
              ENDDO
            ENDDO
            NORMAL_XI=-3
          CASE(GENERATED_MESH_REGULAR_FRONT_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,3)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,1,j),ERR,ERROR)
                SURFACE_NODES(node_counter)=NODE_USER_NUMBER
              ENDDO
            ENDDO
            NORMAL_XI=-2
          CASE(GENERATED_MESH_REGULAR_BACK_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,3)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                NODE_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER(REGULAR_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,SIZE(NIDX,2),j),ERR,ERROR)
                SURFACE_NODES(node_counter)=NODE_USER_NUMBER
              ENDDO
            ENDDO
            NORMAL_XI=2
          CASE DEFAULT
            LOCAL_ERROR="The specified surface type of "//TRIM(NUMBER_TO_VSTRING(SURFACE_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a regular mesh."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Output SURFACE_NODES array is already allocated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Regular mesh object does not have a basis associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Regular mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_REGULAR_SURFACE_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_REGULAR_SURFACE_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_SURFACE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_SURFACE_GET

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_CYLINDER_SURFACE_GET(CYLINDER_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_CYLINDER_TYPE), POINTER :: CYLINDER_MESH !<A pointer to the cylinder mesh object
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GENERATED_MESH_ROUTINES_GeneratedMeshCylinderSurfaces,GENERATED_MESH_ROUTINES
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<On exit, contains the xi direction of the outward pointing normal of the surface
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
    INTEGER(INTG) :: node_counter,i, j, k

    CALL ENTERS("GENERATED_MESH_CYLINDER_SURFACE_GET",ERR,ERROR,*999)

    ! let's go
    IF(ALLOCATED(CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
      NUMBER_OF_ELEMENTS_XI=CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
      IF(ASSOCIATED(CYLINDER_MESH%BASES(MESH_COMPONENT)%PTR)) THEN
        BASIS=>CYLINDER_MESH%BASES(MESH_COMPONENT)%PTR
        IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
          NUMBER_OF_NODES_XI=BASIS%NUMBER_OF_NODES_XIC
          ! build indices first (some of these are dummy arguments)
          CALL GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
              & cylinder_mesh%cylinder_extent,total_number_of_nodes,total_number_of_elements,NIDX,EIDX, &
              & delta,deltai,ERR,ERROR,*999)
          node_counter=0
          SELECT CASE(SURFACE_TYPE)
          CASE(GENERATED_MESH_CYLINDER_INNER_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(1,j,k),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=-1
          CASE(GENERATED_MESH_CYLINDER_OUTER_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,2))*(SIZE(NIDX,3))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO k=1,SIZE(NIDX,3)
              DO j=1,SIZE(NIDX,2)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(SIZE(NIDX,1),j,k),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=1
          CASE(GENERATED_MESH_CYLINDER_TOP_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=3
          CASE(GENERATED_MESH_CYLINDER_BOTTOM_SURFACE)
            ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2))),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
            DO j=1,SIZE(NIDX,2)
              DO i=1,SIZE(NIDX,1)
                node_counter=node_counter+1
                SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(CYLINDER_MESH%GENERATED_MESH,MESH_COMPONENT, &
                    & NIDX(i,j,1),ERR,ERROR)
              ENDDO
            ENDDO
            NORMAL_XI=-3
          CASE DEFAULT
            LOCAL_ERROR="The specified surface type of "//TRIM(NUMBER_TO_VSTRING(SURFACE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Output SURFACE_NODES array is already allocated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Cylinder mesh object does not have a basis associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Cylinder mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_CYLINDER_SURFACE_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CYLINDER_SURFACE_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CYLINDER_SURFACE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_CYLINDER_SURFACE_GET
  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GENERATED_MESH_ELLIPSOID_SURFACE_GET(ELLIPSOID_MESH,MESH_COMPONENT,SURFACE_TYPE,SURFACE_NODES,NORMAL_XI,ERR,ERROR,*)
    ! Argument variables
    TYPE(GENERATED_MESH_ELLIPSOID_TYPE), POINTER :: ELLIPSOID_MESH !<A pointer to the ellipsoid mesh object
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: SURFACE_TYPE !<A constant identifying the type of surface to get \see GENERATED_MESH_ROUTINES_GeneratedMeshEllipsoidSurfaces,GENERATED_MESH_ROUTINES
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: SURFACE_NODES(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: NORMAL_XI !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:),EIDX(:,:,:),CORNER_NODES(:,:,:)
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS_XI(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: NUMBER_OF_NODES_XI(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: total_number_of_nodes,total_number_of_elements
    REAL(DP) :: delta(3),deltai(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: node_counter,i, j, k
    INTEGER(INTG) :: BASIS_COMPONENT

    CALL ENTERS("GENERATED_MESH_ELLIPSOID_SURFACE_GET",ERR,ERROR,*999)

    ! let's go
    IF(ALLOCATED(ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI)) THEN
      NUMBER_OF_ELEMENTS_XI=ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI
      IF(ALLOCATED(ELLIPSOID_MESH%BASES)) THEN

!         !Below, there is an issue:
!         !  BASIS=>ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR does not account for the fact that:
!         !  in 'GENERATED_MESH_ELLIPSOID_CREATE_FINISH' the following is done:
!         !  CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
!         !A temporary work around is the following (although this bug may need to be fixed in several places):
! 
!         IF(MESH_COMPONENT==2) THEN
!           BASIS_COMPONENT = MESH_COMPONENT + 1
!         ELSE
!           BASIS_COMPONENT = MESH_COMPONENT
!         ENDIF
! 
!         IF(ASSOCIATED(ELLIPSOID_MESH%BASES(BASIS_COMPONENT)%PTR)) THEN
!           BASIS=>ELLIPSOID_MESH%BASES(BASIS_COMPONENT)%PTR

        IF(ASSOCIATED(ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR)) THEN
          BASIS=>ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR
          IF(.NOT.ALLOCATED(SURFACE_NODES)) THEN
            NUMBER_OF_NODES_XI=BASIS%NUMBER_OF_NODES_XIC
            ! build indices first (some of these are dummy arguments)
            CALL GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(NUMBER_OF_ELEMENTS_XI,NUMBER_OF_NODES_XI, &
                & ellipsoid_mesh%ellipsoid_extent,total_number_of_nodes,total_number_of_elements,NIDX, &
                & CORNER_NODES,EIDX,delta,deltai,ERR,ERROR,*999)
            node_counter=0

            SELECT CASE(SURFACE_TYPE)

            CASE(GENERATED_MESH_ELLIPSOID_INNER_SURFACE)
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2)-1)+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              j=1
              i=1
              node_counter=node_counter+1
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  NIDX(i,j,1),ERR,ERROR)
              DO j=2,SIZE(NIDX,2)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,j,1)/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,j,1),ERR,ERROR)
                  ELSE
                    node_counter=node_counter-1
                  ENDIF
                ENDDO
              ENDDO
              NORMAL_XI=-3

            CASE(GENERATED_MESH_ELLIPSOID_OUTER_SURFACE)
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,2)-1)+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              j=1
              i=1
              node_counter=node_counter+1
              SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                  & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
              DO j=2,SIZE(NIDX,2)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,j,SIZE(NIDX,3))/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,j,SIZE(NIDX,3)),ERR,ERROR)
                  ELSE
                    node_counter=node_counter-1
                  ENDIF
                ENDDO
              ENDDO
              NORMAL_XI=3

            CASE(GENERATED_MESH_ELLIPSOID_TOP_SURFACE)
              ALLOCATE(SURFACE_NODES((SIZE(NIDX,1))*(SIZE(NIDX,3))),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NODES array.",ERR,ERROR,*999)
              DO k=1,SIZE(NIDX,3)
                DO i=1, SIZE(NIDX,1)
                  node_counter=node_counter+1
                  IF (NIDX(i,SIZE(NIDX,2),k)/=0) THEN
                    SURFACE_NODES(node_counter)=COMPONENT_NODE_TO_USER_NUMBER(ELLIPSOID_MESH%GENERATED_MESH,MESH_COMPONENT, &
                        & NIDX(i,SIZE(NIDX,2),k),ERR,ERROR)
                  ELSE
                    node_counter=node_counter-1
                  ENDIF
                ENDDO
              ENDDO
              NORMAL_XI=2
            CASE DEFAULT
              LOCAL_ERROR="The specified surface type of "//TRIM(NUMBER_TO_VSTRING(SURFACE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Output SURFACE_NODES array is already allocated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Ellipsoid mesh object does not have the first basis associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Ellipsoid mesh object does not have bases allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Ellipsoid mesh object does not have number of elements property specified.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("GENERATED_MESH_ELLIPSOID_SURFACE_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ELLIPSOID_SURFACE_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ELLIPSOID_SURFACE_GET")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_SURFACE_GET

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given regular mesh (Not to be called by user)
  SUBROUTINE GENERATED_MESH_REGULAR_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XIC,MAXIMUM_EXTENT, &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAi,ERR,ERROR,*)
    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XIC(3) !<Number of nodes per element in each xi direction for this basis
    REAL(DP),INTENT(IN) :: MAXIMUM_EXTENT(3)         !<width, length and height of regular mesh
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_NODES    !<On exit, contains total number of nodes in regular mesh component
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS !<On exit, contains total number of elements in regular mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: NIDX(:,:,:)  !<Mapping array to find a node number for a given (x,y,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: EIDX(:,:,:)  !<Mapping array to find an element number for a given (x,y,z)
    REAL(DP),INTENT(OUT) :: DELTA(3)  !<Step sizes in each of (x,y,z) for elements
    REAL(DP),INTENT(OUT) :: DELTAi(3) !<Step sizes in each of (x,y,z) for node (identical to DELTA if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) !<Total number of nodes per element in each xi direction for this basis

    CALL ENTERS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES",ERR,ERROR,*999)

    IF(.NOT.ALLOCATED(NIDX)) THEN
      IF(.NOT.ALLOCATED(EIDX)) THEN
        ! calculate DELTA and DELTAi
        DELTA(1)=MAXIMUM_EXTENT(1)/NUMBER_ELEMENTS_XI(1)
        DELTA(2)=MAXIMUM_EXTENT(2)/NUMBER_ELEMENTS_XI(2)
        DELTA(3)=MAXIMUM_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
        DO xi_idx=1,3
          DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
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

    CALL EXITS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES")
    RETURN
999 CALL ERRORS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_BUILD_NODE_INDICES")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_BUILD_NODE_INDICES

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given cylinder (Not to be called by user)
  SUBROUTINE GENERATED_MESH_CYLINDER_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XIC,CYLINDER_EXTENT, &
    & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,EIDX,DELTA,DELTAi,ERR,ERROR,*)
    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XIC(3) !<Number of nodes per element in each xi direction (basis property)
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
          DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XIC(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XIC(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(2)=TOTAL_NUMBER_NODES_XI(2)-1 ! theta loops around so slightly different
        !TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

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

  !
  !================================================================================================================================
  !

  !>Calculate the mesh topology information for a given ellipsoid (Not to be called by user)
  SUBROUTINE GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES(NUMBER_ELEMENTS_XI,NUMBER_OF_NODES_XI,ELLIPSOID_EXTENT, &
    & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NIDX,CORNER_NODES,EIDX,DELTA,DELTAi,ERR,ERROR,*)
    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: NUMBER_ELEMENTS_XI(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: NUMBER_OF_NODES_XI(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: ELLIPSOID_EXTENT(4)         !< long axis, short axis, wall thickness, top angle
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_NODES    !<On exit, contains total number of nodes in ellipsoid mesh
    INTEGER(INTG),INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS !<On exit, contains total number of elements in ellipsoid mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: NIDX(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: CORNER_NODES(:,:,:) ! Returns the array of corner nodes numbered
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: EIDX(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: DELTA(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: DELTAi(3) !<Step sizes in each of (r,theta,z) for node (identical to DELTA if 2 nodes per xi direction)
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string

    ! Local variables
    INTEGER(INTG) :: xi_idx,ne1,ne2,ne3,nn1,nn2,nn3,tn1,tn2,tn3,NN,NE
    INTEGER(INTG) :: TOTAL_NUMBER_NODES_XI(3) ! total number of nodes in each xi direction

    CALL ENTERS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES",ERR,ERROR,*999)

    ! Can skip most of the testing as this subroutine is only to be called by
    ! GENERATED_MESH_ELLIPSOID_CREATE_FINISH, which tests the input params.
    IF(.NOT.ALLOCATED(NIDX)) THEN
      IF(.NOT.ALLOCATED(EIDX)) THEN
        ! calculate DELTA and DELTAi
        DELTA(1)=TWOPI/NUMBER_ELEMENTS_XI(1)
        DELTA(2)=(PI-ELLIPSOID_EXTENT(4))/NUMBER_ELEMENTS_XI(2)
        DELTA(3)=ELLIPSOID_EXTENT(3)/NUMBER_ELEMENTS_XI(3)
        DO xi_idx=1,3
          DELTAi(xi_idx)=DELTA(xi_idx)/(NUMBER_OF_NODES_XI(xi_idx)-1)
        ENDDO

        ! calculate total elements and nodes
        DO xi_idx=1,3
          TOTAL_NUMBER_NODES_XI(xi_idx)=(NUMBER_OF_NODES_XI(xi_idx)-1)*NUMBER_ELEMENTS_XI(xi_idx)+1
        ENDDO
        TOTAL_NUMBER_NODES_XI(1)=TOTAL_NUMBER_NODES_XI(1)-1 ! circumferential loops around so slightly different
        TOTAL_NUMBER_OF_ELEMENTS=PRODUCT(NUMBER_ELEMENTS_XI)

        ! calculate NIDX first
        ALLOCATE(NIDX(TOTAL_NUMBER_NODES_XI(1),TOTAL_NUMBER_NODES_XI(2),TOTAL_NUMBER_NODES_XI(3)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NIDX array.",ERR,ERROR,*999)
        ALLOCATE(CORNER_NODES(NUMBER_ELEMENTS_XI(1),NUMBER_ELEMENTS_XI(2)+1,NUMBER_ELEMENTS_XI(3)+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NIDX array.",ERR,ERROR,*999)

        !nn: node number inside element in certain direction
        !ne: element number in certain direction
        !tn: global node number in certain direction
        !NN: Node counter
        !Due to one more corner node than elements in transmural direction, first shell is taken separatly
        NN=0
        ne3=1
        nn3=1
        !Due to one more corner node than elements in longitudinal direction, apex elements are taken separatly
        ne2=1
        nn2=1
        ne1=1
        nn1=1
        !apex nodes
        NN=NN+1
        tn3=1
        tn2=1
        tn1=1
        NIDX(tn1,tn2,tn3)=NN
        CORNER_NODES(ne1,ne2,ne3)=NN
        DO ne2=1,NUMBER_ELEMENTS_XI(2)
          DO nn2=2,(NUMBER_OF_NODES_XI(2))
            tn2=tn2+1
            tn1=0
            DO ne1=1,NUMBER_ELEMENTS_XI(1)
              DO nn1=1,(NUMBER_OF_NODES_XI(1)-1)
                tn1=tn1+1
                NN=NN+1
                NIDX(tn1,tn2,tn3)=NN
                IF ((nn1==1).AND.(nn2==NUMBER_OF_NODES_XI(2))) THEN
                  CORNER_NODES(ne1,ne2+1,ne3)=NN
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO ne3=1,NUMBER_ELEMENTS_XI(3)
          DO nn3=2,NUMBER_OF_NODES_XI(3)
            ne2=1
            nn2=1
            ne1=1
            nn1=1
            !apex nodes
            NN=NN+1
            tn3=tn3+1
            tn2=1
            tn1=1
            NIDX(tn1,tn2,tn3)=NN
            IF (nn3==NUMBER_OF_NODES_XI(3)) THEN
              CORNER_NODES(ne1,ne2,ne3+1)=NN
            ENDIF
            DO ne2=1,NUMBER_ELEMENTS_XI(2)
              DO nn2=2,(NUMBER_OF_NODES_XI(2))
                tn2=tn2+1
                tn1=0
                DO ne1=1,NUMBER_ELEMENTS_XI(1)
                  DO nn1=1,(NUMBER_OF_NODES_XI(1)-1)
                    tn1=tn1+1
                    NN=NN+1
                    NIDX(tn1,tn2,tn3)=NN
                    IF ((nn1==1).AND.(nn3==NUMBER_OF_NODES_XI(3)).AND.(nn2==NUMBER_OF_NODES_XI(2))) THEN
                      CORNER_NODES(ne1,ne2+1,ne3+1)=NN
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
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

    CALL EXITS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ELLIPSOID_BUILD_NODE_INDICES

  !
  !================================================================================================================================
  !

  !>Calculates the user node numbers for an array of nodes numbered using one basis
  SUBROUTINE COMPONENT_NODES_TO_USER_NUMBERS(GENERATED_MESH,BASIS_INDEX,NODE_COMPONENT_NUMBERS, &
      & NODE_USER_NUMBERS,ERR,ERROR,*)
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH   !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_COMPONENT_NUMBERS(:)  !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: NODE_USER_NUMBERS(:)    !<On return, the corresponding user numbers
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !local variables
    INTEGER(INTG) :: node_idx

    CALL ENTERS("COMPONENT_NODES_TO_USER_NUMBERS",ERR,ERROR,*999)

    IF(SIZE(NODE_USER_NUMBERS)==SIZE(NODE_COMPONENT_NUMBERS)) THEN
      DO node_idx=1,SIZE(NODE_COMPONENT_NUMBERS)
        NODE_USER_NUMBERS(node_idx)=COMPONENT_NODE_TO_USER_NUMBER(GENERATED_MESH,BASIS_INDEX, &
            NODE_COMPONENT_NUMBERS(node_idx),ERR,ERROR)
      ENDDO
    ELSE
      CALL FLAG_ERROR("NODE_COMPONENT_NUMBERS and NODE_USER_NUMBERS arrays have different sizes.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN
999 CALL ERRORS("COMPONENT_NODES_TO_USER_NUMBERS",ERR,ERROR)
    CALL EXITS("COMPONENT_NODES_TO_USER_NUMBERS")
    RETURN 1
  END SUBROUTINE COMPONENT_NODES_TO_USER_NUMBERS

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis
  FUNCTION COMPONENT_NODE_TO_USER_NUMBER(GENERATED_MESH,BASIS_INDEX,NODE_COMPONENT_NUMBER,ERR,ERROR)
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH  !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                  !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_COMPONENT_NUMBER           !<The node number for this component basis
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !function variable
    INTEGER(INTG) :: COMPONENT_NODE_TO_USER_NUMBER !<On return, the corresponding user node number

    !local variables
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,basis_idx,ni,REMAINDER,REMAINDER2,TEMP_TERM,NUM_CORNER_NODES,NODE_OFFSET,BASIS_NUM_NODES
    INTEGER(INTG) :: POS(3),POS2(3),CORNER_NODE_FACTOR(3),BASIS_NODE_FACTOR(3),BASIS_ELEMENT_FACTOR(3),NUM_PREVIOUS_CORNERS,STEP
    INTEGER(INTG), POINTER :: NUMBER_OF_ELEMENTS_XI(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    LOGICAL :: CORNER_NODE,FINISHED_COUNT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("COMPONENT_NODE_TO_USER_NUMBER",ERR,ERROR,*999)

    NULLIFY(BASIS)
    NULLIFY(BASES)
    NUM_CORNER_NODES=1
    REMAINDER=NODE_COMPONENT_NUMBER-1 !use zero based numbering
    REMAINDER2=NODE_COMPONENT_NUMBER-1
    COMPONENT_NODE_TO_USER_NUMBER=0
    POS=0
    POS2=0

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%REGULAR_MESH%BASES)
          NUM_DIMS=GENERATED_MESH%REGULAR_MESH%MESH_DIMENSION
          BASES=>GENERATED_MESH%REGULAR_MESH%BASES
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
        CALL FLAG_ERROR("The regular mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%CYLINDER_MESH%BASES)
          NUM_DIMS=GENERATED_MESH%CYLINDER_MESH%MESH_DIMENSION
          BASES=>GENERATED_MESH%CYLINDER_MESH%BASES
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
        CALL FLAG_ERROR("The cylinder mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%ELLIPSOID_MESH)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%ELLIPSOID_MESH%BASES)
          NUM_DIMS=GENERATED_MESH%ELLIPSOID_MESH%MESH_DIMENSION
          BASES=>GENERATED_MESH%ELLIPSOID_MESH%BASES
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%ELLIPSOID_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
        CALL FLAG_ERROR("The ellipsoid mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh generated type of "// &
            & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      IF(BASIS_INDEX<=NUM_BASES) THEN
        IF(NUM_BASES==1) THEN
          !If is the only basis, don't do anything
          COMPONENT_NODE_TO_USER_NUMBER=NODE_COMPONENT_NUMBER
        ELSE
          TEMP_TERM=1
          NUM_CORNER_NODES=1
          DO ni=1,NUM_DIMS
            NUM_CORNER_NODES=NUM_CORNER_NODES*(NUMBER_OF_ELEMENTS_XI(ni)+1)
            CORNER_NODE_FACTOR(ni)=1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*(NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              CORNER_NODE_FACTOR(ni)=CORNER_NODE_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            CORNER_NODE_FACTOR(3)=CORNER_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)-1
            NUM_CORNER_NODES=NUM_CORNER_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_ELEMENTS_XI(3)+1)
          ELSE IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            CORNER_NODE_FACTOR(3)=CORNER_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)-NUMBER_OF_ELEMENTS_XI(2)
            CORNER_NODE_FACTOR(2)=CORNER_NODE_FACTOR(2)-1
            NUM_CORNER_NODES=NUM_CORNER_NODES-(NUMBER_OF_ELEMENTS_XI(2)+1)*(NUMBER_OF_ELEMENTS_XI(3)+1)- &
                & (NUMBER_OF_ELEMENTS_XI(1)-1)*(NUMBER_OF_ELEMENTS_XI(3)+1)
          ENDIF
          NODE_OFFSET=NUM_CORNER_NODES
          IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            !Every second mesh component is the collapsed node version
            STEP=2
          ELSE
            STEP=1
          ENDIF
          DO basis_idx=1,BASIS_INDEX-1,STEP
            BASIS=>BASES(basis_idx)%PTR
            BASIS_NUM_NODES=1
            DO ni=1,NUM_DIMS
              BASIS_NUM_NODES=BASIS_NUM_NODES*(NUMBER_OF_ELEMENTS_XI(ni)*(BASIS%NUMBER_OF_NODES_XIC(ni)-1)+1)
            ENDDO
            !Adjust for other mesh types
            IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
              BASIS_NUM_NODES=BASIS_NUM_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(BASIS%NUMBER_OF_nodes_xic(1)-1)* &
                  & (NUMBER_OF_ELEMENTS_XI(3)+1)*(BASIS%NUMBER_OF_nodes_xic(3)-1)
            ELSE IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
              BASIS_NUM_NODES=BASIS_NUM_NODES-(NUMBER_OF_ELEMENTS_XI(2)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)+1)* &
                & (NUMBER_OF_ELEMENTS_XI(3)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)+1)- &
                & (NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1)* &
                & (NUMBER_OF_ELEMENTS_XI(3)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)+1)
            ENDIF
            NODE_OFFSET=NODE_OFFSET+BASIS_NUM_NODES-NUM_CORNER_NODES
          ENDDO
          BASIS=>BASES(BASIS_INDEX)%PTR
          TEMP_TERM=1
          DO ni=1,NUM_DIMS
            BASIS_NODE_FACTOR(ni)=1
            BASIS_ELEMENT_FACTOR(ni)=BASIS%NUMBER_OF_NODES_XIC(ni)-1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*((BASIS%NUMBER_OF_NODES_XIC(ni-1)-1)*NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              BASIS_NODE_FACTOR(ni)=BASIS_NODE_FACTOR(ni)*TEMP_TERM
              BASIS_ELEMENT_FACTOR(ni)=BASIS_ELEMENT_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            !subtract nodes along line where y wraps around
            BASIS_NODE_FACTOR(3)=BASIS_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(1)-1)+1)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)
          ELSE IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            !subtract missing nodes at apex
            BASIS_NODE_FACTOR(3)=BASIS_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)+1
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(1)-1)+1)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)
            !subtract nodes along line where x wraps around
            BASIS_NODE_FACTOR(3)=BASIS_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(2)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)-1
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(2)*(BASIS%NUMBER_OF_NODES_XIC(2)-1)-1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(3)-1)
            BASIS_NODE_FACTOR(2)=BASIS_NODE_FACTOR(2)-1
            BASIS_ELEMENT_FACTOR(2)=BASIS_ELEMENT_FACTOR(2)-(BASIS%NUMBER_OF_NODES_XIC(2)-1)
          ENDIF
          !Work out if we have a corner node, otherwise add node numbers used by corners and
          !previous basis interpolations and subtract number of corner nodes used before the
          !given component node number to get the user number
          CORNER_NODE=.TRUE.
          IF(NUM_DIMS>2) THEN
            POS(3)=REMAINDER/BASIS_NODE_FACTOR(3)
            POS2(3)=REMAINDER2/BASIS_ELEMENT_FACTOR(3)
            REMAINDER=MOD(REMAINDER,BASIS_NODE_FACTOR(3))
            REMAINDER2=MOD(REMAINDER2,BASIS_ELEMENT_FACTOR(3))
            IF(MOD(POS(3),BASIS%NUMBER_OF_NODES_XIC(3)-1)/=0) CORNER_NODE=.FALSE.
          ENDIF
          IF(NUM_DIMS>1) THEN
            IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
              !Need to account for missing nodes at apex
              IF(REMAINDER>0) THEN
                REMAINDER=REMAINDER+NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1
                REMAINDER2=REMAINDER2+NUMBER_OF_ELEMENTS_XI(1)*(BASIS%NUMBER_OF_NODES_XIC(1)-1)-1
              ENDIF
            ENDIF
            POS(2)=REMAINDER/BASIS_NODE_FACTOR(2)
            POS2(2)=REMAINDER2/BASIS_ELEMENT_FACTOR(2)
            REMAINDER=MOD(REMAINDER,BASIS_NODE_FACTOR(2))
            REMAINDER2=MOD(REMAINDER2,BASIS_ELEMENT_FACTOR(2))
            IF(MOD(POS(2),BASIS%NUMBER_OF_NODES_XIC(2)-1)/=0) CORNER_NODE=.FALSE.
          ENDIF
          POS(1)=REMAINDER/BASIS_NODE_FACTOR(1)
          POS2(1)=REMAINDER2/BASIS_ELEMENT_FACTOR(1)
          IF(MOD(POS(1),BASIS%NUMBER_OF_NODES_XIC(1)-1)/=0) CORNER_NODE=.FALSE.
          IF(CORNER_NODE) THEN
            COMPONENT_NODE_TO_USER_NUMBER=POS2(1)*CORNER_NODE_FACTOR(1)+POS2(2)*CORNER_NODE_FACTOR(2)+ &
                & POS2(3)*CORNER_NODE_FACTOR(3)
            IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE.AND.POS2(2)/=0) THEN
              !Subtract off non-existent nodes at apex
              COMPONENT_NODE_TO_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER-(NUMBER_OF_ELEMENTS_XI(1)-1)
            ENDIF
            COMPONENT_NODE_TO_USER_NUMBER=COMPONENT_NODE_TO_USER_NUMBER+1
          ELSE
            !subtract previous corner nodes from node offset
            NUM_PREVIOUS_CORNERS=0
            FINISHED_COUNT=.FALSE.
            IF(NUM_DIMS>2) THEN
              IF(MOD(POS(3),BASIS%NUMBER_OF_NODES_XIC(3)-1)/=0) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*(POS2(3)+1)
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*POS2(3)
              ENDIF
            ENDIF
            IF(NUM_DIMS>1 .AND. FINISHED_COUNT.NEQV..TRUE.) THEN
              IF(MOD(POS(2),BASIS%NUMBER_OF_NODES_XIC(2)-1)/=0) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*(POS2(2)+1)
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*POS2(2)
              ENDIF
              IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS-(NUMBER_OF_ELEMENTS_XI(1)-1)
              ENDIF
            ENDIF
            IF(FINISHED_COUNT.NEQV..TRUE.) THEN
              IF(MOD(POS(1),BASIS%NUMBER_OF_NODES_XIC(1)-1)/=0) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(1)*(POS2(1)+1)
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(1)*POS2(1)
              ENDIF
            ENDIF
            NODE_OFFSET=NODE_OFFSET-NUM_PREVIOUS_CORNERS
            COMPONENT_NODE_TO_USER_NUMBER=NODE_OFFSET+NODE_COMPONENT_NUMBER
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="Mesh component must be less than or equal to "//(NUMBER_TO_VSTRING(NUM_BASES,"*",ERR,ERROR))// &
            & " but it is "//(NUMBER_TO_VSTRING(BASIS_INDEX,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("COMPONENT_NODE_TO_USER_NUMBER")
    RETURN
999 CALL ERRORS("COMPONENT_NODE_TO_USER_NUMBER",ERR,ERROR)
    CALL EXITS("COMPONENT_NODE_TO_USER_NUMBER")
    RETURN
  END FUNCTION COMPONENT_NODE_TO_USER_NUMBER

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis
  FUNCTION USER_NUMBER_TO_COMPONENT_NODE(GENERATED_MESH,BASIS_INDEX,NODE_USER_NUMBER,ERR,ERROR)
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH        !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: BASIS_INDEX                     !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: NODE_USER_NUMBER                !<The corresponding user node number
    INTEGER(INTG) :: ERR          !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !function variable
    INTEGER(INTG) :: USER_NUMBER_TO_COMPONENT_NODE !<On return, the node number for this component basis
    !local variables
    INTEGER(INTG) :: NUM_BASES,NUM_DIMS,basis_idx,ni,REMAINDER,TEMP_TERM,NUM_CORNER_NODES,NODE_OFFSET,BASIS_NUM_NODES
    INTEGER(INTG) :: POS(3),CORNER_NODE_FACTOR(3),BASIS_ELEMENT_FACTOR(3),NUM_PREVIOUS_CORNERS
    INTEGER(INTG), POINTER :: NUMBER_OF_ELEMENTS_XI(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:)
    LOGICAL :: FINISHED_COUNT,OFF_EDGE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("USER_NUMBER_TO_COMPONENT_NODE",ERR,ERROR,*999)

    NULLIFY(BASIS)
    NULLIFY(BASES)
    NUM_CORNER_NODES=1
    REMAINDER=NODE_USER_NUMBER-1 !use zero based numbering
    POS=0

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%REGULAR_MESH%BASES)
          NUM_DIMS=GENERATED_MESH%REGULAR_MESH%MESH_DIMENSION
          BASES=>GENERATED_MESH%REGULAR_MESH%BASES
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
        CALL FLAG_ERROR("The regular mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%CYLINDER_MESH)) THEN
          NUM_BASES=SIZE(GENERATED_MESH%CYLINDER_MESH%BASES)
          NUM_DIMS=GENERATED_MESH%CYLINDER_MESH%MESH_DIMENSION
          BASES=>GENERATED_MESH%CYLINDER_MESH%BASES
          NUMBER_OF_ELEMENTS_XI=>GENERATED_MESH%CYLINDER_MESH%NUMBER_OF_ELEMENTS_XI
        ELSE
          CALL FLAG_ERROR("The cylinder mesh for this generated mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        !Ellipsoid mesh does things a bit differently and this function isn't needed
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The generated mesh generated type of "// &
            & TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%GENERATED_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      IF(BASIS_INDEX<=NUM_BASES) THEN
        IF(NUM_BASES==1) THEN
          !If is the only basis, don't do anything
          USER_NUMBER_TO_COMPONENT_NODE=NODE_USER_NUMBER
        ELSE
          TEMP_TERM=1
          NUM_CORNER_NODES=1
          DO ni=1,NUM_DIMS
            NUM_CORNER_NODES=NUM_CORNER_NODES*(NUMBER_OF_ELEMENTS_XI(ni)+1)
            CORNER_NODE_FACTOR(ni)=1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*(NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              CORNER_NODE_FACTOR(ni)=CORNER_NODE_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            CORNER_NODE_FACTOR(3)=CORNER_NODE_FACTOR(3)-NUMBER_OF_ELEMENTS_XI(1)-1
            NUM_CORNER_NODES=NUM_CORNER_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(NUMBER_OF_ELEMENTS_XI(3)+1)
          ENDIF
          NODE_OFFSET=NUM_CORNER_NODES
          DO basis_idx=1,BASIS_INDEX-1
            BASIS=>BASES(basis_idx)%PTR
            BASIS_NUM_NODES=1
            DO ni=1,NUM_DIMS
              BASIS_NUM_NODES=BASIS_NUM_NODES*(NUMBER_OF_ELEMENTS_XI(ni)*(BASIS%NUMBER_OF_NODES_XIC(ni)-1)+1)
            ENDDO
            !Adjust for other mesh types
            IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
              BASIS_NUM_NODES=BASIS_NUM_NODES-(NUMBER_OF_ELEMENTS_XI(1)+1)*(BASIS%NUMBER_OF_nodes_xic(1)-1)* &
                  & (NUMBER_OF_ELEMENTS_XI(3)+1)*(BASIS%NUMBER_OF_nodes_xic(3)-1)
            ENDIF
            NODE_OFFSET=NODE_OFFSET+BASIS_NUM_NODES-NUM_CORNER_NODES
          ENDDO
          BASIS=>BASES(BASIS_INDEX)%PTR
          TEMP_TERM=1
          DO ni=1,NUM_DIMS
            BASIS_ELEMENT_FACTOR(ni)=BASIS%NUMBER_OF_NODES_XIC(ni)-1
            IF(ni>1) THEN
              TEMP_TERM=TEMP_TERM*((BASIS%NUMBER_OF_NODES_XIC(ni-1)-1)*NUMBER_OF_ELEMENTS_XI(ni-1)+1)
              BASIS_ELEMENT_FACTOR(ni)=BASIS_ELEMENT_FACTOR(ni)*TEMP_TERM
            ENDIF
          ENDDO
          !Adjust for other mesh types
          IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
            !subtract nodes along line where y wraps around
            BASIS_ELEMENT_FACTOR(3)=BASIS_ELEMENT_FACTOR(3)-(NUMBER_OF_ELEMENTS_XI(1)* &
                & (BASIS%NUMBER_OF_NODES_XIC(1)-1)+1)*(BASIS%NUMBER_OF_NODES_XIC(3)-1)
          ENDIF
          IF(NODE_USER_NUMBER<=NUM_CORNER_NODES) THEN
            !we have a node on a corner
            IF(NUM_DIMS>2) THEN
              POS(3)=REMAINDER/CORNER_NODE_FACTOR(3)
              REMAINDER=MOD(REMAINDER,CORNER_NODE_FACTOR(3))
            ENDIF
            IF(NUM_DIMS>1) THEN
              POS(2)=REMAINDER/CORNER_NODE_FACTOR(2)
              REMAINDER=MOD(REMAINDER,CORNER_NODE_FACTOR(2))
            ENDIF
            POS(1)=REMAINDER/CORNER_NODE_FACTOR(1)
            USER_NUMBER_TO_COMPONENT_NODE=POS(1)*BASIS_ELEMENT_FACTOR(1)+POS(2)*BASIS_ELEMENT_FACTOR(2)+ &
                & POS(3)*BASIS_ELEMENT_FACTOR(3)
            USER_NUMBER_TO_COMPONENT_NODE=USER_NUMBER_TO_COMPONENT_NODE+1
          ELSE IF(NODE_USER_NUMBER>NODE_OFFSET) THEN
            REMAINDER=REMAINDER-NODE_OFFSET
            DO ni=1,NUM_DIMS
              BASIS_ELEMENT_FACTOR(ni)=BASIS_ELEMENT_FACTOR(ni)-CORNER_NODE_FACTOR(ni)
            ENDDO
            NUM_PREVIOUS_CORNERS=0
            FINISHED_COUNT=.FALSE.
            OFF_EDGE=.FALSE.
            IF(NUM_DIMS>2) THEN
              IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_CYLINDER_MESH_TYPE.AND. &
                  & (MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3)) > BASIS_ELEMENT_FACTOR(2)*NUMBER_OF_ELEMENTS_XI(2)-1)) THEN
                OFF_EDGE=.TRUE.
              ELSE IF(GENERATED_MESH%GENERATED_TYPE==GENERATED_MESH_REGULAR_MESH_TYPE.AND. &
                  & MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3)) > (BASIS_ELEMENT_FACTOR(2)*NUMBER_OF_ELEMENTS_XI(2)+ &
                  & BASIS_ELEMENT_FACTOR(1)*NUMBER_OF_ELEMENTS_XI(1)-1)) THEN
                OFF_EDGE=.TRUE.
              ENDIF
              IF(OFF_EDGE) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*(1+REMAINDER/BASIS_ELEMENT_FACTOR(3))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3))
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(3)*(REMAINDER/BASIS_ELEMENT_FACTOR(3))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(3))
              ENDIF
            ENDIF
            IF(NUM_DIMS>1 .AND. FINISHED_COUNT.NEQV..TRUE.) THEN
              IF(MOD(REMAINDER,BASIS_ELEMENT_FACTOR(2)) > &
                  & BASIS_ELEMENT_FACTOR(1)*NUMBER_OF_ELEMENTS_XI(1)-1) THEN
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*(1+REMAINDER/BASIS_ELEMENT_FACTOR(2))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(2))
                FINISHED_COUNT=.TRUE.
              ELSE
                NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(2)*(REMAINDER/BASIS_ELEMENT_FACTOR(2))
                REMAINDER=MOD(REMAINDER,BASIS_ELEMENT_FACTOR(2))
              ENDIF
            ENDIF
            IF(FINISHED_COUNT.NEQV..TRUE.) THEN
              NUM_PREVIOUS_CORNERS=NUM_PREVIOUS_CORNERS+CORNER_NODE_FACTOR(1)*(REMAINDER/BASIS_ELEMENT_FACTOR(1))+1
            ENDIF
            NODE_OFFSET=NODE_OFFSET-NUM_PREVIOUS_CORNERS
            USER_NUMBER_TO_COMPONENT_NODE=NODE_USER_NUMBER-NODE_OFFSET
          ELSE
            CALL FLAG_ERROR("Invalid node number specified.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="Mesh component must be less than or equal to "//(NUMBER_TO_VSTRING(NUM_BASES,"*",ERR,ERROR))// &
            & " but it is "//(NUMBER_TO_VSTRING(BASIS_INDEX,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("USER_NUMBER_TO_COMPONENT_NODE")
    RETURN
999 CALL ERRORS("USER_NUMBER_TO_COMPONENT_NODE",ERR,ERROR)
    CALL EXITS("USER_NUMBER_TO_COMPONENT_NODE")
    RETURN
  END FUNCTION USER_NUMBER_TO_COMPONENT_NODE

  !
  !================================================================================================================================
  !

END MODULE GENERATED_MESH_ROUTINES

