!> \file
!> $Id: generated_mesh_routines.f90 9 2007-05-15 13:52:02Z cpb $
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

!> This module handles all generated mesh routines.
MODULE GENERATED_MESH_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE COORDINATE_ROUTINES
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshTypes GENERATED_MESH_ROUTINES::GeneratedMeshTypes 
  !> \brief Generated mesh types.
  !> \see GENERATED_MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_MESH_TYPE=1 !<A regular generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_POLAR_MESH_TYPE=2 !<A polar generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_FRACTAL_TREE_MESH_TYPE=3 !<A fractal tree generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables
  TYPE(GENERATED_MESHES_TYPE), TARGET :: GENERATED_MESHES

  !Interfaces
  
  INTERFACE GENERATED_MESH_BASIS_SET
    MODULE PROCEDURE GENERATED_MESH_BASIS_SET_NUMBER
    MODULE PROCEDURE GENERATED_MESH_BASIS_SET_PTR
  END INTERFACE !GENERATED_MESH_BASIS_SET
  
  INTERFACE GENERATED_MESH_DESTROY
    MODULE PROCEDURE GENERATED_MESH_DESTROY_NUMBER
    MODULE PROCEDURE GENERATED_MESH_DESTROY_PTR
  END INTERFACE !GENERATED_MESH_DESTROY
  
  INTERFACE GENERATED_MESH_EXTENT_SET
    MODULE PROCEDURE GENERATED_MESH_EXTENT_SET_NUMBER
    MODULE PROCEDURE GENERATED_MESH_EXTENT_SET_PTR
  END INTERFACE !GENERATED_MESH_EXTENT_SET
  
  INTERFACE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET
    MODULE PROCEDURE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER
    MODULE PROCEDURE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR
  END INTERFACE !GENERATED_MESH_NUMBER_OF_ELEMENTS_SET
  
  INTERFACE GENERATED_MESH_ORIGIN_SET
    MODULE PROCEDURE GENERATED_MESH_ORIGIN_SET_NUMBER
    MODULE PROCEDURE GENERATED_MESH_ORIGIN_SET_PTR
  END INTERFACE !GENERATED_MESH_ORIGIN_SET
  
  INTERFACE GENERATED_MESH_TYPE_SET
    MODULE PROCEDURE GENERATED_MESH_TYPE_SET_NUMBER
    MODULE PROCEDURE GENERATED_MESH_TYPE_SET_PTR
  END INTERFACE !GENERATED_MESH_TYPE_SET
  
  
  PUBLIC GENERATED_MESHES_INITIALISE,GENERATED_MESHES_FINALISE,GENERATED_MESH_CREATE_START,GENERATED_MESH_CREATE_FINISH 
  
  PUBLIC GENERATED_MESH_BASIS_SET,GENERATED_MESH_EXTENT_SET,GENERATED_MESH_NUMBER_OF_ELEMENTS_SET,GENERATED_MESH_ORIGIN_SET, &
    & GENERATED_MESH_TYPE_SET
    
  PUBLIC GENERATED_MESH_BASIS_GET,GENERATED_MESH_EXTENT_GET,GENERATED_MESH_NUMBER_OF_ELEMENTS_GET,GENERATED_MESH_ORIGIN_GET,&
    & GENERATED_MESH_TYPE_GET


CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Gets the basis of a generated mesh.
  FUNCTION GENERATED_MESH_BASIS_GET(GENERATED_MESH,ERR,ERROR)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the basis of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function result
    TYPE(BASIS_TYPE) :: GENERATED_MESH_BASIS_GET !<The basis of mesh to generate \see GENERATED_MESH_ROUTINES_GeneratedMeshBasis,GENERATED_MESH_ROUTINES
    !Local Variables

    CALL ENTERS("GENERATED_MESH_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
          GENERATED_MESH_BASIS_GET=GENERATED_MESH%REGULAR_MESH%BASIS
        ELSE
          CALL FLAG_ERROR("Regular generated mesh is not associated",ERR,ERROR,*999)
        END IF
      CASE DEFAULT
        CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_BASIS_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_BASIS_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_BASIS_GET")
    RETURN
  END FUNCTION GENERATED_MESH_BASIS_GET

  !
  !================================================================================================================================
  ! 
  
  !>Set the basis of the generated mesh. 
  SUBROUTINE GENERATED_MESH_BASIS_SET_NUMBER(USER_NUMBER,BASIS,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_BASIS_SET_NUMBER
    !###  Description:
    !###    Sets/changes the basis of a generated mesh identified by a USER_NUMBER.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH

    CALL ENTERS("GENERATED_MESH_BASIS_SET_NUMBER",ERR,ERROR,*999)

    CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
    CALL GENERATED_MESH_BASIS_SET_PTR(GENERATED_MESH,BASIS,ERR,ERROR,*999)
    
    CALL EXITS("GENERATED_MESH_BASIS_SET_NUMBER")
    RETURN
999 CALL ERRORS("GENERATED_MESH_BASIS_SET_NUMBER",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_BASIS_SET_NUMBER")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_BASIS_SET_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the basis of a generated mesh.
  SUBROUTINE GENERATED_MESH_BASIS_SET_PTR(GENERATED_MESH,BASIS,ERR,ERROR,*)

     !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the basis of
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The basis of mesh to generate \see GENERATED_MESH_ROUTINES_GeneratedMeshBasis,GENERATED_MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_BASIS_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BASIS)) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
            IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
              GENERATED_MESH%REGULAR_MESH%BASIS=>BASIS 
            ELSE
              CALL FLAG_ERROR("Regular generated mesh is not associated",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
          END SELECT
          
        ELSE
          CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_BASIS_SET_PTR")
    RETURN
999 CALL ERRORS("GENERATED_MESH_BASIS_SET_PTR",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_BASIS_SET_PTR")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_BASIS_SET_PTR
  
  !
  !================================================================================================================================
  !
  
  !>Finishes the creation of a generated mesh.
  SUBROUTINE GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,MESH,MESH_USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to finish the creation of
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the generated mesh to finish the creation of
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<The mesh's user number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR,*999)
    
    NULLIFY(MESH)
    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%REGION)) THEN
        IF(ASSOCIATED(GENERATED_MESH%REGION%COORDINATE_SYSTEM)) THEN
          SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE) 
            CALL GENERATED_MESH_REGULAR_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*999)
            !CALL MESH_CREATE_START(MESH_USER_NUMBER,GENERATED_MESH%REGION,GENERATED_MESH%REGULAR_MESH%MESH_DIMENSION,MESH,ERR,ERROR,*999)
            MESH=>GENERATED_MESH%MESH
            MESH%GENERATED_MESH=>GENERATED_MESH
          CASE DEFAULT
            CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
          END SELECT
          MESH=>GENERATED_MESH%MESH
        ELSE
          CALL FLAG_ERROR("Region coordinate system is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
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

  !>Starts the creation of a generated mesh.
  SUBROUTINE GENERATED_MESH_CREATE_START(USER_NUMBER,REGION,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the generated mesh to create
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the generated mesh on
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: generated_mesh_idx
    TYPE(GENERATED_MESH_TYPE), POINTER :: NEW_GENERATED_MESH
    TYPE(GENERATED_MESH_PTR_TYPE), POINTER :: NEW_GENERATED_MESHES(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    NULLIFY(NEW_GENERATED_MESH)
    NULLIFY(NEW_GENERATED_MESHES)

    CALL ENTERS("GENERATED_MESH_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(GENERATED_MESH)) THEN
      LOCAL_ERROR="Generated Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
        & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(GENERATED_MESH%REGION%USER_NUMBER,"*",ERR,ERROR))
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ELSE
      IF(ASSOCIATED(REGION)) THEN
        ALLOCATE(NEW_GENERATED_MESH,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new generated mesh.",ERR,ERROR,*999)
        ! Initialise generated mesh
        CALL GENERATED_MESH_INITIALISE(NEW_GENERATED_MESH,ERR,ERROR,*999)
        !Set default generated mesh values
        NEW_GENERATED_MESH%USER_NUMBER=USER_NUMBER
        NEW_GENERATED_MESH%GLOBAL_NUMBER=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES+1
        NEW_GENERATED_MESH%REGION=>REGION
        
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
        GENERATED_MESH=>NEW_GENERATED_MESH
	  ELSE
        CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
      ENDIF
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_START")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_START",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_START")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Destroys a generated mesh.
  SUBROUTINE GENERATED_MESH_DESTROY_NUMBER(USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN):: USER_NUMBER !<The user number of generated mesh to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: generated_mesh_position
    LOGICAL :: FOUND
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GENERATED_MESH_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESHES%GENERATED_MESHES)) THEN
      !Find the generated mesh identified by the user number
      FOUND=.FALSE.
      generated_mesh_position=0
      DO WHILE(generated_mesh_position<GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES.AND..NOT.FOUND)
        generated_mesh_position=generated_mesh_position+1
        IF(GENERATED_MESHES%GENERATED_MESHES(generated_mesh_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
      ENDDO
      
      IF(FOUND) THEN
        GENERATED_MESH=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_position)%PTR
        !Destroy all the generated mesh components
        CALL GENERATED_MESH_DESTROY(GENERATED_MESH,ERR,ERROR,*999)
      ELSE
        LOCAL_ERROR="Generated mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))//" has not been created."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated meshes is not associated.",ERR,ERROR,*999)
    ENDIF    
 
    CALL EXITS("GENERATED_MESH_DESTROY_NUMBER")
    RETURN
999 CALL ERRORS("GENERATED_MESH_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_DESTROY_NUMBER")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh.
  SUBROUTINE GENERATED_MESH_DESTROY_PTR(GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: generated_mesh_idx,generated_mesh_position
    TYPE(GENERATED_MESH_PTR_TYPE), POINTER :: NEW_GENERATED_MESHES(:)

    CALL ENTERS("GENERATED_MESH_DESTROY_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
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
              GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR%GLOBAL_NUMBER=GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR%GLOBAL_NUMBER-1
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
        CALL FLAG_ERROR("Generated meshes are not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    END IF
 
    CALL EXITS("GENERATED_MESH_DESTROY_PTR")
    RETURN
998 IF(ASSOCIATED(NEW_GENERATED_MESHES)) DEALLOCATE(NEW_GENERATED_MESHES)
999 CALL ERRORS("GENERATED_MESH_DESTROY_PTR",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_DESTROY_PTR")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_DESTROY_PTR
  
  !
  !================================================================================================================================
  !

  !>Gets the extent of a generated mesh.
  FUNCTION GENERATED_MESH_EXTENT_GET(GENERATED_MESH,ERR,ERROR)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the type of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Funtion result
    REAL(DP), POINTER :: GENERATED_MESH_EXTENT_GET(:)
    !Local Variables

    CALL ENTERS("GENERATED_MESH_EXTENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        ALLOCATE(GENERATED_MESH_EXTENT_GET(SIZE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate generated mesh extent.",ERR,ERROR,*999)
        GENERATED_MESH_EXTENT_GET=GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT
      CASE DEFAULT
        CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_EXTENT_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_EXTENT_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_EXTENT_GET")
    RETURN   
  END FUNCTION GENERATED_MESH_EXTENT_GET
  
  !
  !================================================================================================================================
  ! 
   
  !>Sets/changes the max extent of a generated mesh. 
  SUBROUTINE GENERATED_MESH_EXTENT_SET_NUMBER(USER_NUMBER,MAX_EXTENT,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_EXTENT_SET_NUMBER
    !###  Description:
    !###    Sets/changes the max extent of a generated mesh identified by a USER_NUMBER.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    REAL(DP), INTENT(IN) :: MAX_EXTENT(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH

    CALL ENTERS("GENERATED_MESH_EXTENT_SET_NUMBER",ERR,ERROR,*999)

    CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
    CALL GENERATED_MESH_EXTENT_SET_PTR(GENERATED_MESH,MAX_EXTENT,ERR,ERROR,*999)
    
    CALL EXITS("GENERATED_MESH_EXTENT_SET_NUMBER")
    RETURN
999 CALL ERRORS("GENERATED_MESH_EXTENT_SET_NUMBER",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_EXTENT_SET_NUMBER")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_EXTENT_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the extent of a generated mesh.
  SUBROUTINE GENERATED_MESH_EXTENT_SET_PTR(GENERATED_MESH,MAX_EXTENT,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: MAX_EXTENT(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_EXTENT_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
            ALLOCATE(GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT(SIZE(MAX_EXTENT)),STAT=ERR)
            GENERATED_MESH%REGULAR_MESH%MAXIMUM_EXTENT=MAX_EXTENT 
          ELSE
            CALL FLAG_ERROR("Regular generated mesh is not associated",ERR,ERROR,*999)
          END IF
        CASE DEFAULT
          CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_EXTENT_SET_PTR")
    RETURN
999 CALL ERRORS("GENERATED_MESH_EXTENT_SET_PTR",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_EXTENT_SET_PTR")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_EXTENT_SET_PTR

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
      GENERATED_MESH%USER_NUMBER=0
      GENERATED_MESH%GLOBAL_NUMBER=0
      NULLIFY(GENERATED_MESH%REGION)
      GENERATED_MESH%GENERATED_TYPE=0
      NULLIFY(GENERATED_MESH%REGULAR_MESH)
      NULLIFY(GENERATED_MESH%MESH)
      !Default to a regular mesh.
      CALL GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
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

  !>Gets the extent of a generated mesh.
  FUNCTION GENERATED_MESH_NUMBER_OF_ELEMENTS_GET(GENERATED_MESH,ERR,ERROR)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function result
    INTEGER(INTG), POINTER :: GENERATED_MESH_NUMBER_OF_ELEMENTS_GET(:)
    !Local Variables

    CALL ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        ALLOCATE(GENERATED_MESH_NUMBER_OF_ELEMENTS_GET(SIZE(GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate generated mesh number of elements.",ERR,ERROR,*999)
        GENERATED_MESH_NUMBER_OF_ELEMENTS_GET=GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI
      CASE DEFAULT
        CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN 
  END FUNCTION GENERATED_MESH_NUMBER_OF_ELEMENTS_GET
  
  !
  !================================================================================================================================
  ! 
    
  SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER(USER_NUMBER,NUMBER_OF_ELEMENTS_XI,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_BASIS_SET_NUMBER
    !###  Description:
    !###    Sets/changes the basis of a generated mesh identified by a USER_NUMBER.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH

    CALL ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER",ERR,ERROR,*999)

    CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
    CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR(GENERATED_MESH,NUMBER_OF_ELEMENTS_XI,ERR,ERROR,*999)
    
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER")
    RETURN
999 CALL ERRORS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_NUMBER
                 

  !
  !================================================================================================================================
  !

  !>Sets/changes the extent of a generated mesh.
  SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR(GENERATED_MESH,NUMBER_OF_ELEMENTS_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
        CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
          IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN 
            ALLOCATE(GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(SIZE(NUMBER_OF_ELEMENTS_XI)),STAT=ERR)
            GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI=NUMBER_OF_ELEMENTS_XI 
          ELSE
            CALL FLAG_ERROR("Regular generated mesh is not associated",ERR,ERROR,*999)
          END IF
        CASE DEFAULT
          CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR")
    RETURN
999 CALL ERRORS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_NUMBER_OF_ELEMENTS_SET_PTR
  
  !
  !================================================================================================================================
  !

  !>Get the origin of a generated mesh.
  FUNCTION GENERATED_MESH_ORIGIN_GET(GENERATED_MESH,ERR,ERROR)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to get the type of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function result
    REAL(DP), POINTER :: GENERATED_MESH_ORIGIN_GET(:)
    !Local Variables

    CALL ENTERS("GENERATED_MESH_ORIGIN_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      SELECT CASE(GENERATED_MESH%GENERATED_TYPE)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        ALLOCATE(GENERATED_MESH_ORIGIN_GET(SIZE(GENERATED_MESH%REGULAR_MESH%ORIGIN)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate generated mesh origin.",ERR,ERROR,*999)
        GENERATED_MESH_ORIGIN_GET=GENERATED_MESH%REGULAR_MESH%ORIGIN 
      CASE DEFAULT
        CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_ORIGIN_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ORIGIN_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ORIGIN_GET")
    RETURN
  END FUNCTION GENERATED_MESH_ORIGIN_GET
  
  !
  !================================================================================================================================
  ! 
  
  !>Sets/changes the origin of a generated mesh.  
  SUBROUTINE GENERATED_MESH_ORIGIN_SET_NUMBER(USER_NUMBER,ORIGIN,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_ORIGIN_SET_NUMBER
    !###  Description:
    !###    Sets/changes the basis of a generated mesh identified by a USER_NUMBER.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    REAL(DP), INTENT(IN) :: ORIGIN(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH

    CALL ENTERS("GENERATED_MESH_ORIGIN_SET_NUMBER",ERR,ERROR,*999)

    CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
    CALL GENERATED_MESH_ORIGIN_SET_PTR(GENERATED_MESH,ORIGIN,ERR,ERROR,*999)
    
    CALL EXITS("GENERATED_MESH_ORIGIN_SET_NUMBER")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ORIGIN_SET_NUMBER",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ORIGIN_SET_NUMBER")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_ORIGIN_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a generated mesh.
  SUBROUTINE GENERATED_MESH_ORIGIN_SET_PTR(GENERATED_MESH,ORIGIN,ERR,ERROR,*)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: ORIGIN(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESH_ORIGIN_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished",ERR,ERROR,*999)
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
            CALL FLAG_ERROR("Regular generated mesh is not associated",ERR,ERROR,*999)
          END IF
        CASE DEFAULT
          CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_ORIGIN_SET_PTR")
    RETURN
999 CALL ERRORS("GENERATED_MESH_ORIGIN_SET_PTR",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_ORIGIN_SET_PTR")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_ORIGIN_SET_PTR
  
  !
  !================================================================================================================================
  ! 
  
  !>Start to create the regular generated mesh type
  SUBROUTINE GENERATED_MESH_REGULAR_CREATE_FINISH(GENERATED_MESH,MESH_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: origin_idx
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), ALLOCATABLE :: NUMBER_ELEMENTS_XI(:)
    TYPE(REGION_TYPE), POINTER :: REGION 
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG) :: ni,ne,ne1,ne2,ne3,NN,nn1,nn2,nn3,np,TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_ELEMENTS_XI(3), &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    
    CALL ENTERS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
        REGULAR_MESH=>GENERATED_MESH%REGULAR_MESH
        
        ! Validating
        IF(ASSOCIATED(GENERATED_MESH%REGION)) THEN
          REGION=>GENERATED_MESH%REGION
          IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
            !TODO is regular type only for COORDINATE_RECTANGULAR_CARTESIAN_TYPE? 
            !If that, should we use IF rather than select?
            SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              !Determine the coordinate system and create the regular mesh for that system
              REGULAR_MESH%MESH_DIMENSION=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              NUMBER_OF_DIMENSIONS=REGULAR_MESH%MESH_DIMENSION
              IF (.NOT.ALLOCATED(REGULAR_MESH%ORIGIN)) THEN
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
                        IF(.NOT.ALL(NUMBER_ELEMENTS_XI>0)) CALL FLAG_ERROR("Must have 1 or more elements in all directions",ERR,ERROR,*999)
                        IF(.NOT.ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED)) CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented",ERR,ERROR,*999)
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
                        CALL NODES_CREATE_START(TOTAL_NUMBER_OF_NODES,REGION,NODES,ERR,ERROR,*999)
	                    !Create the mesh
	                    CALL MESH_CREATE_START(MESH_USER_NUMBER,REGION,SIZE(NUMBER_ELEMENTS_XI,1),GENERATED_MESH%MESH,ERR,ERROR,*999)
	                    !Create the elements
	                    CALL MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH%MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
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
	                                    ne=ne+(ne3-1)*TOTAL_NUMBER_ELEMENTS_XI(1)*TOTAl_NUMBER_ELEMENTS_XI(2)
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
	                    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(GENERATED_MESH%MESH,1,ERR,ERROR,*999)
	                    !Finish the mesh
	                    CALL MESH_CREATE_FINISH(REGION,GENERATED_MESH%MESH,ERR,ERROR,*999)

                      ELSE
                        CALL FLAG_ERROR("The number of xi directions of the given basis does not match the size of &
                          &the number of elements for the mesh",ERR,ERROR,*999)
                      ENDIF  
                    CASE(BASIS_SIMPLEX_TYPE)                  
                      CALL FLAG_ERROR("Regular meshes with simplex basis types is not implemented",ERR,ERROR,*999)
		            CASE DEFAULT
		              CALL FLAG_ERROR("Basis type is either invalid or not implemented",ERR,ERROR,*999)
		            END SELECT
                  ELSE
                    CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                    &the maximum extent",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The number of dimensions of the given regular mesh does not match the size of &
                  &the origin",ERR,ERROR,*999)
              ENDIF
            CASE default
              CALL FLAG_ERROR("Coordinate type is either invalid or not implemented",ERR,ERROR,*999)
            END SELECT
		  ELSE
		    CALL FLAG_ERROR("Coordiate System is not associated",ERR,ERROR,*999)
		  ENDIF
		ELSE
          CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Regular mesh is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN
    ! TODO invalidate other associations
999 CALL ERRORS("GENERATED_MESH_REGULAR_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_CREATE_FINISH
  
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
    
    CALL ENTERS("GENERATED_MESH_REGULAR_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(ASSOCIATED(GENERATED_MESH%REGULAR_MESH)) THEN
        CALL FLAG_ERROR("regular mesh type is already associated for this generated mesh",ERR,ERROR,*999)
      ELSE
        ALLOCATE(GENERATED_MESH%REGULAR_MESH,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated regular generated mesh",ERR,ERROR,*999)
        GENERATED_MESH%REGULAR_MESH%GENERATED_MESH=>GENERATED_MESH
        GENERATED_MESH%GENERATED_TYPE=GENERATED_MESH_REGULAR_MESH_TYPE
        NULLIFY(GENERATED_MESH%REGULAR_MESH%BASIS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN
    ! TODO invalidate other associations
999 CALL ERRORS("GENERATED_MESH_REGULAR_INITIALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_REGULAR_INITIALISE")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_REGULAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the type of a generated mesh.
  FUNCTION GENERATED_MESH_TYPE_GET(GENERATED_MESH,ERR,ERROR)

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function result
    INTEGER(INTG) :: GENERATED_MESH_TYPE_GET !<The type of mesh to generate \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
    !Local Variables

    CALL ENTERS("GENERATED_MESH_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      GENERATED_MESH_TYPE_GET=GENERATED_MESH%GENERATED_TYPE
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_TYPE_GET")
    RETURN
999 CALL ERRORS("GENERATED_MESH_TYPE_GET",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_TYPE_GET")
    RETURN 
  END FUNCTION GENERATED_MESH_TYPE_GET
  
  !
  !================================================================================================================================
  ! 
   
  !>Sets/changes the type of a generated mesh.  
  SUBROUTINE GENERATED_MESH_TYPE_SET_NUMBER(USER_NUMBER,GENERATED_TYPE,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_TYPE_SET_NUMBER
    !###  Description:
    !###    Sets/changes the type of a generated mesh identified by a USER_NUMBER.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    INTEGER(INTG), INTENT(IN) :: GENERATED_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH

    CALL ENTERS("GENERATED_MESH_TYPE_SET_NUMBER",ERR,ERROR,*999)

    CALL GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,GENERATED_MESH,ERR,ERROR,*999)
    CALL GENERATED_MESH_TYPE_SET_PTR(GENERATED_MESH,GENERATED_TYPE,ERR,ERROR,*999)
    
    CALL EXITS("GENERATED_MESH_TYPE_SET_NUMBER")
    RETURN
999 CALL ERRORS("GENERATED_MESH_TYPE_SET_NUMBER",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_TYPE_SET_NUMBER")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh.
  SUBROUTINE GENERATED_MESH_TYPE_SET_PTR(GENERATED_MESH,GENERATED_TYPE,ERR,ERROR,*)

     !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: GENERATED_TYPE !<The type of mesh to generate \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     INTEGER(INTG) :: OLD_GENERATED_TYPE

    CALL ENTERS("GENERATED_MESH_TYPE_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      IF(GENERATED_MESH%GENERATED_MESH_FINISHED) THEN
        CALL FLAG_ERROR("Generated mesh has been finished",ERR,ERROR,*999)
      ELSE
        OLD_GENERATED_TYPE=GENERATED_MESH%GENERATED_TYPE
        IF(OLD_GENERATED_TYPE/=GENERATED_TYPE) THEN
	      SELECT CASE(GENERATED_TYPE)
	      CASE(GENERATED_MESH_REGULAR_MESH_TYPE) 
	        CALL GENERATED_MESH_REGULAR_INITIALISE(GENERATED_MESH,ERR,ERROR,*999)
	      CASE DEFAULT
	        CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
	      END SELECT
	      
	      SELECT CASE(OLD_GENERATED_TYPE)
          CASE(GENERATED_MESH_REGULAR_MESH_TYPE) 
            CALL GENERATED_MESH_REGULAR_FINALISE(GENERATED_MESH%REGULAR_MESH,ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Generated mesh type is either invalid or not implemented",ERR,ERROR,*999)
          END SELECT
	    ENDIF
	  ENDIF
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_TYPE_SET_PTR")
    RETURN
999 CALL ERRORS("GENERATED_MESH_TYPE_SET_PTR",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_TYPE_SET_PTR")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_TYPE_SET_PTR
  
  !
  !================================================================================================================================
  !
  
  !>Finds and returns in generated mesh a pointer to that identified by USER_NUMBER. If no generated mesh with that number exits left nullified.
  SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND(USER_NUMBER,GENERATED_MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: generated_mesh_idx

    CALL ENTERS("GENERATED_MESH_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      CALL FLAG_ERROR("Generated mesh is already associated.",ERR,ERROR,*999)
    ELSE
      generated_mesh_idx=1
      DO WHILE(generated_mesh_idx<=GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES.AND..NOT.ASSOCIATED(GENERATED_MESH))
        IF(GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          GENERATED_MESH=>GENERATED_MESHES%GENERATED_MESHES(generated_mesh_idx)%PTR
        ELSE
          generated_mesh_idx=generated_mesh_idx+1
        ENDIF
      ENDDO
    ENDIF
    
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("MESH_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE GENERATED_MESH_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises all generated meshes and deallocates all memory.
  SUBROUTINE GENERATED_MESHES_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::USER_NUMBER

    CALL ENTERS("GENERATED_MESHES_FINALISE",ERR,ERROR,*999)

    DO WHILE(GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES>0)
      USER_NUMBER=GENERATED_MESHES%GENERATED_MESHES(1)%PTR%USER_NUMBER
      CALL GENERATED_MESH_DESTROY(USER_NUMBER,ERR,ERROR,*999)
    ENDDO !generated_mesh_idx
    
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
  SUBROUTINE GENERATED_MESHES_INITIALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GENERATED_MESHES_INITIALISE",ERR,ERROR,*999)

    GENERATED_MESHES%NUMBER_OF_GENERATED_MESHES=0
    NULLIFY(GENERATED_MESHES%GENERATED_MESHES)
    
    CALL EXITS("GENERATED_MESHES_INITIALISE")
    RETURN
999 CALL ERRORS("GENERATED_MESHES_INITIALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESHES_INITIALISE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESHES_INITIALISE
  
  !
  !================================================================================================================================
  ! 

END MODULE GENERATED_MESH_ROUTINES
